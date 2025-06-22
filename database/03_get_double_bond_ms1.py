import re 
import json


import numpy as np 
import pandas as pd 
from rdkit import Chem
from rdkit.Chem import Descriptors
from pathos import multiprocessing as mp
import matchms

from tqdm import tqdm 


H_MASS = 1.008
O_MASS = 15.9994

def simple_parallel(
    input_list, function, n_threads=12, timeout=4000, max_retries=3
):
    cpus = min(mp.cpu_count(), n_threads)
    pool = mp.ProcessingPool(nodes=cpus)
    results = list(tqdm(pool.imap(function, input_list), total=len(input_list)))
    return results


def chunked_parallel(
    input_list, function, chunck_size = 10000, n_threads=12, timeout=4000, max_retries=3
):
    def batch_func(list_inputs):
        outputs = []
        for i in list_inputs:
            outputs.append(function(i))
        return outputs
    
    chunked_list = [
        input_list[i : i + chunck_size] for i in range(0, len(input_list), chunck_size)
    ]
    list_outputs = simple_parallel(
        chunked_list,
        batch_func,
        n_threads=n_threads,
        timeout=timeout,
        max_retries=max_retries,
    )
    # Unroll
    full_output = [j for i in list_outputs for j in i]
    return full_output


def normlize_smiles(smi: str):
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)


def smiles_to_MolWt(smi: str):
    mol = Chem.MolFromSmiles(smi)
    try:
        mol_weight = Descriptors.MolWt(mol)
    except:
        mol_weight = 0.0
    return mol_weight
        

def get_ozonolysis_delta_mass(metadata: dict):
    
    DB_PATTERN = r'\d+:\d+\([\d,ZE]+\)'
    chain_shorthand = metadata['chain_shorthand']
    lipid = metadata["db_shorthand"]
    matches = re.findall(DB_PATTERN, lipid)
    if len(matches) == 0:
        return {'chain_shorthand': chain_shorthand, 'lipid_name': lipid, 'delta_mass':[]}
    
    mz_shift_list = []

    for s in matches:
        s = re.sub('[ZE]', '', s)
        num_fa, num_db, pos_db_delta = re.split(r'[:\(]', s)
        pos_db_delta = re.sub(r'\)', '', pos_db_delta)
        pos_db_delta = re.split(r'\,', pos_db_delta)
        pos_db_omega  = [int(num_fa) - int(db) for db in pos_db_delta]
        pos_db_omega.sort(reverse=False)
        for db in pos_db_omega:
            smiles = "C" * db
            aldehyde_mz_shift = -1.0 * smiles_to_MolWt(smiles) + 2 * H_MASS * (pos_db_omega.index(db) + 1) + O_MASS
            criegee_mz_shift = aldehyde_mz_shift + O_MASS
            mz_shift_list.extend([aldehyde_mz_shift, criegee_mz_shift])
        
    ozonide_mz_shift = 3 * O_MASS
    mz_shift_list.append(ozonide_mz_shift)
    
    mz_shift_list = list(set(mz_shift_list))
    mz_shift_list.sort(reverse=False)
    
    return {'chain_shorthand': chain_shorthand, 'lipid_name': lipid, 'delta_mass': mz_shift_list}
        

if __name__ == "__main__":
    in_file = "clean_data/lipidmaps_chain_shorthand.tsv"
    out_file = "clean_data/ozonolysis_delta_mass.json"

    metadata_list = pd.read_csv(in_file, sep='\t').to_dict("records")
    delta_mass_list = chunked_parallel(metadata_list, get_ozonolysis_delta_mass, chunck_size=1000)
    
    with open(out_file, 'w') as f:
        json.dump(delta_mass_list, f)
import re 
import json


import numpy as np 
import pandas as pd
import anndata as ad
from matchms.importing import load_from_json, load_from_mgf, load_from_msp


def parse_mzmine_peak_table(file: str) -> ad.AnnData:
    """
        Parse peak table from mzmine.
    Args:
    
    Return:
    
    Raises:    
    
    """
    df = pd.read_csv(file)
    
    feature_cols = [
        'id', 'rt', 'mz', 
        'spectral_db_matches:compound_name', 'spectral_db_matches:ion_adduct', 'spectral_db_matches:cosine_score',
        'spectral_db_matches:mol_formula', 'spectral_db_matches:smiles', 'spectral_db_matches:inchi'
    ]  
    data_cols = df.columns[df.columns.str.match('datafile.*area')]
    samples = [re.findall(r':(.*?):', col)[0] for col in data_cols]
    
    X = df[data_cols].fillna(0).to_numpy().T
    vars = df[feature_cols]
    vars = vars.sort_values(['rt', 'mz'], ascending=[True, True])
    vars.columns = vars.columns.str.replace(r'spectral_db_matches:', '')
    vars = vars.fillna('')
    vars = vars.set_index('id')
    vars.index = vars.index.astype(str)
    #vars.index = pd.Series(['F' + str(i).rjust(4, '0') for i in range(len(vars))])
    obs = pd.DataFrame(data={'sample': samples})
    #obs.index=pd.Series(['S' + str(i).rjust(3, '0') for i in range(len(obs))])
    obs.index = obs.index.astype(str)

    return ad.AnnData(X, var=vars, obs=obs)


def load_ms1_peak_table(filename, mode = "mzmine"):
    match mode:
        case "mzmine":
            df_ms1 = load_mzmine_peak_table(filename)
        case _:
            pass
    return df_ms1
    

def load_mzmine_peak_table(filename):
    df = pd.read_csv(filename)
    feature_cols = [
        'id', 'rt', 'mz', 
        'spectral_db_matches:compound_name', 'spectral_db_matches:ion_adduct', 'spectral_db_matches:cosine_score',
        'spectral_db_matches:mol_formula', 'spectral_db_matches:smiles', #'spectral_db_matches:inchi'
    ] 
    data_cols = df.columns[df.columns.str.match('datafile.*area')].to_list()
    df_data = df[data_cols].fillna(0)

    df_vars = df[feature_cols]
    df_vars = df_vars.sort_values(['rt', 'mz'], ascending=[True, True])
    df_vars.columns = df_vars.columns.str.replace(r'spectral_db_matches:', '')
    df_vars = df_vars.rename(columns={'ion_adduct': 'adduct'})
    df_vars = df_vars.fillna('')
    
    df_ms1 = pd.concat([df_vars, df_data],axis=1)
    return df_ms1


def parse_ms1_peaktable(file: str):
    pass

def parse_ms2_spectra(file: str):
    pass


def parse_msp_file():
    pass 

def parse_mgf_file():
    pass

def parse_json_file():
    pass


def load_spec_file(filename, metadata_harmonization=True):
    spectrum_list = list(load_from_mgf(filename, metadata_harmonization))
    spectra = {spectrum.metadata["feature_id"] : spectrum for spectrum in spectrum_list}
    return spectra


def load_structure_database(filename):
    """"""
    df = pd.read_csv(filename, sep='\t')
    df = df[['compound_name', 'adduct', 'chain_shorthand', 'db_shorthand']]
    return df 


def load_ozonolysis_prodcut_database(filename):
    with open(filename, 'r') as f:
        delta_mass_list = json.load(f)
    df = pd.DataFrame().from_records(delta_mass_list)
    df = df[df['delta_mass'].apply(lambda x: x != [] if isinstance(x, list) else False)]
    return df 


def load_reference_database(structure_file, ozid_file):
    df_structure = load_structure_database(structure_file)
    df_ozid = load_ozonolysis_prodcut_database(ozid_file)
    df_reference = pd.merge(df_structure, df_ozid, how="inner", on=["chain_shorthand"])
    return df_reference
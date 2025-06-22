import pandas as pd 
import matchms
from matchms.importing import load_from_msp, load_from_mgf

from pathos import multiprocessing as mp
from tqdm import tqdm


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
    # Adding it here fixes somessetting disrupted elsewhere
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


def extract_spectrum_metadata(spectrum: matchms.Spectrum):
    metadata = spectrum.metadata 
    spec_metadata = {
        "spectrum_id": metadata.get("spectrum_id", ""),
        "compound_name": metadata.get("compound_name", ""),
        "smiles": metadata.get("smiles", ""),
        "inchikey": metadata.get("inchikey", ""),
        "precursor_mz": metadata.get("precursor_mz", 0.0),
        "adduct": metadata.get("adduct", ""),
        "formula": metadata.get("formula", ""),
        "ionmode": metadata.get("ionmode", ""),
        "instrument": metadata.get("instrument", ""),
        "instrument_type": metadata.get("instrument_type", ""),
        "collision_energy": metadata.get("collision_energy", ""),
    }
    return spec_metadata


if __name__ == "__main__":
    lipid_pos_file = "raw_data/MSDIAL-TandemMassSpectralAtlas-VS69-Lipid-Pos.msp"
    lipid_neg_file = "raw_data/MSDIAL-TandemMassSpectralAtlas-VS69--Lipid-Neg.msp"
    out_file_neg = "clean_data/ms-dial_ms2_metadata_neg.tsv"
    out_file_pos = "clean_data/ms-dial_ms2_metadata_pos.tsv"

    print("Loading Spectra!")
    lipid_pos_spectra = list(load_from_msp(lipid_pos_file))
    lipid_neg_spectra = list(load_from_msp(lipid_neg_file))
    print("Done!")

    records_pos_metadata = chunked_parallel(lipid_pos_spectra, extract_spectrum_metadata, chunck_size=1000)
    df_pos_metadata = pd.DataFrame.from_records(records_pos_metadata)
    df_pos_metadata.to_csv(out_file_pos, sep='\t', index_label="ID")
    

    records_neg_metadata = chunked_parallel(lipid_neg_spectra, extract_spectrum_metadata, chunck_size=1000)
    df_neg_metadata = pd.DataFrame.from_records(records_neg_metadata)
    df_neg_metadata.to_csv(out_file_neg, sep='\t', index_label="ID")

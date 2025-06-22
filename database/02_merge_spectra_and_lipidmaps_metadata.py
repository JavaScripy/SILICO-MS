import re
import pandas as pd


SP_PATTERN = r'\b(Cer|AHexCer|HexCer|SL|SM)\b'
OH_ORIGIN_PATTERN_LIST = [r";O", r";2O", r";3O"]
OH_REPLACE_PATTERN_LIST = ["m", "d", "t"]
OH_PATTERN_LIST = list(zip(OH_ORIGIN_PATTERN_LIST, OH_REPLACE_PATTERN_LIST))


def normlize_sp_chain_shorthand(shorthand):
    for origin, replacement in OH_PATTERN_LIST:
        hydroxyl_pattern_list = re.findall(origin, shorthand)
        if len(hydroxyl_pattern_list) > 0:
            shorthand = re.sub(origin + "/", "_", shorthand)
            shorthand = re.sub(" ", "(" + replacement, shorthand)
            shorthand = shorthand + ")"
            break
    return shorthand

def normlize_normal_chain_shorthand(shorthand):
    shorthand = re.sub(" ", "(", shorthand)
    shorthand = shorthand + ")"
    return shorthand


def get_chain_shorthand(row):
    chain_shorthand  = row["compound_name"]

    sp_pattern_list = re.findall(SP_PATTERN, chain_shorthand)
    
    if len(sp_pattern_list) > 0:
        chain_shorthand = normlize_sp_chain_shorthand(chain_shorthand)
    else:
        chain_shorthand = normlize_normal_chain_shorthand(chain_shorthand)
    return chain_shorthand


def merge_metadata(structure_file, spctra_file, out_file):
    df_structure = pd.read_csv(structure_file, sep='\t')
    df_spectra = pd.read_csv(spctra_file, sep='\t')
    df_spectra  = df_spectra.drop(["smiles", "inchikey", "formula"], axis=1)
    df_spectra ["chain_shorthand"] = df_spectra .apply(get_chain_shorthand, axis=1)

    df_toatal = pd.merge(df_spectra, df_structure, how="inner", on=["chain_shorthand"])
    df_toatal.to_csv(out_file, sep='\t', index=None)

if __name__ == "__main__":
    structure_file = "clean_data/lipidmaps_chain_shorthand.tsv"
    # neg mode 
    spctra_file = "clean_data/ms-dial_ms2_metadata_neg.tsv"
    out_file = "clean_data/reference_ms2_neg.tsv"

    merge_metadata(structure_file, spctra_file, out_file)
    
    # pos mode
    spctra_file = "clean_data/ms-dial_ms2_metadata_pos.tsv"
    out_file = "clean_data/reference_ms2_pos.tsv"
    
    merge_metadata(structure_file, spctra_file, out_file)
    
    

    


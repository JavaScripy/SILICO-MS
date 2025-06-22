import re
import numpy as np 
import pandas as pd


FA_PATTEN =  r'\d+:\d+'
DB_PATTERN = r"\([\d,ZE]+\)"
SUBGROUP_PATTERN = r'\([\dMeOH]+\)'


def normlize_shorthand(shorthand):
    shorthand = shorthand.strip(" ")
    shorthand = shorthand.replace(r"[R]", "")
    shorthand = shorthand.replace(r"[rac]", "")
    shorthand = shorthand.replace(r"[iso2]", "")
    shorthand = shorthand.replace(r"[iso3]", "")
    shorthand = shorthand.replace(r"[iso6]", "")

    shorthand = shorthand.replace(r"/0:0", "")
    shorthand = shorthand.replace(r"(0:0/", "")
    shorthand = shorthand.replace(r"/", "_")
    # CL
    shorthand = shorthand.replace(r"1'-[", "")
    shorthand = shorthand.replace(r",3'-[", "_")
    shorthand = shorthand.replace(r"]", "")
    return shorthand


def normlize_chain_shorthand(shorthand):
    shorthand = shorthand.strip(" ")
    shorthand = re.sub(DB_PATTERN, "", shorthand) 
    shorthand = re.sub(SUBGROUP_PATTERN, "", shorthand)
    return shorthand


def normlize_db_shorthand(shorthand):
    shorthand = shorthand.strip(" ")
    shorthand = re.sub(r"(\d)[ZE]", r"\1", shorthand) 
    return shorthand


def get_shorthand(row):
    common_name = row["common name"]
    if not isinstance(common_name, str):
        shorthand = row["species shorthand"]
    else:
        fa_pattern_list = re.findall(FA_PATTEN, common_name)
        if len(fa_pattern_list) > 0:
            shorthand = normlize_shorthand(common_name)
        else:
            shorthand = row["species shorthand"]
    
    return shorthand


def get_chain_shorthand(row):
    chain_shorthand = get_shorthand(row)
    if not isinstance(chain_shorthand, str):
        chain_shorthand = row["species shorthand"]
    else:
        db_pattern_list = re.findall(DB_PATTERN, chain_shorthand)
        subgroup_pattern_list = re.findall(SUBGROUP_PATTERN, chain_shorthand)
        if len(db_pattern_list) + len(subgroup_pattern_list) > 0:
            chain_shorthand = normlize_chain_shorthand(chain_shorthand)
        else:
            chain_shorthand = row["shorthand"]
    return chain_shorthand


def get_db_shorthand(row):
    db_shorthand = get_shorthand(row)
    if not isinstance(db_shorthand, str):
        db_shorthand = row["species shorthand"]
    else:
        db_pattern_list = re.findall(DB_PATTERN, db_shorthand)
        if len(db_pattern_list) > 0:
            db_shorthand = normlize_db_shorthand(db_shorthand)
        else:
            db_shorthand = row["shorthand"]
    return db_shorthand
            


def get_merged_df(id_file, shorthand_file):
    df_id = pd.read_csv(id_file, sep='\t')
    df_shorthand = pd.read_csv(shorthand_file, sep='\t')
    df = pd.merge(df_id, df_shorthand, how="inner", on=['lm_id'])
    df = df.drop(['obsolete_id', 'systematic name'], axis=1)
    return df 

if __name__ == "__main__":
    id_file = "raw_data/lipidmaps_ids_cc0.tsv"
    shorthand_file = "raw_data/lipidmaps_shorthand.tsv"
    out_file = "clean_data/lipidmaps_chain_shorthand.tsv"

    df = get_merged_df(id_file, shorthand_file)
    
    df["shorthand"] = df.apply(get_shorthand, axis=1)
    df["chain_shorthand"] = df.apply(get_chain_shorthand, axis=1)
    df["db_shorthand"] = df.apply(get_db_shorthand, axis=1)
    df = df.drop(["common name"], axis=1)
    
    df.to_csv(out_file, sep='\t', index=None)
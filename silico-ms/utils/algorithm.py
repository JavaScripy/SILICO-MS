from typing import List, Dict


import pandas as pd
import networkx as nx 
from matchms import Spectrum 
from matchms.similarity import CosineGreedy, ModifiedCosine, NeutralLossesCosine, CosineHungarian


def get_similarity_type(score_type: str = "None",
                        mz_tol: float = 0.1):
    match score_type:
        case "NIST-LC":
            # reference:
            similarity_score = CosineGreedy(tolerance=mz_tol, mz_power=1.3, intensity_power=0.53)
        case "NIST-GC":
            # reference:
            similarity_score = CosineGreedy(tolerance=mz_tol, mz_power=3.0, intensity_power=0.6)
        case "SQRT":
            # reference:
            similarity_score = CosineGreedy(tolerance=mz_tol, mz_power=0.0, intensity_power=0.5)
        case "MassBank":
            # reference:
            similarity_score = CosineGreedy(tolerance=mz_tol, mz_power=2.0, intensity_power=0.5)
        case "None":
            # reference:
            similarity_score = CosineGreedy(tolerance=mz_tol, mz_power=0.0, intensity_power=1.0)
        case "ModifiedCosine":
            # reference
            similarity_score = ModifiedCosine(tolerance=mz_tol, mz_power=1.3, intensity_power=0.53)
        case "NeutralLossesCosine":
            # reference
            similarity_score = NeutralLossesCosine(tolerance=mz_tol, mz_power=1.3, intensity_power=0.53)
        case "CosineHungarian":
            # reference
            similarity_score = CosineHungarian(tolerance=mz_tol, mz_power=1.3, intensity_power=0.53)
        case _:
            print("Not matched method!")
    
    return similarity_score


def get_rt_similarity(feature1: dict, 
                      feature2: dict, 
                      tolerance: float = 0.5, 
                      mode: str = "absolute"):
    """edited from MetDNA.
    Args:
        feature1: query
        feature2: reference
    """
    rt1 = feature1["rt"]
    rt2 = feature2["rt"]
    
    match mode:
        case "relative":
            rt_diff = abs(rt1 - rt2) * 100 / rt2
        case "absolute":
            rt_diff = abs(rt1 - rt2)
        case _:
            pass
    
    if rt_diff < tolerance:
        rt_score = 1 - abs(rt1 - rt2) / tolerance
    else:
        rt_score = 0.0
    return rt_score


def get_precursor_mz_similarity(feature1: dict, 
                                feature2: dict, 
                                tolerance = 5, 
                                mode = "ppm"):
    """edited from MetDNA.
    Args:
        feature1: query
        feature2: reference
    """
    mz1 = feature1["mz"]
    mz2 = feature2["mz"]   
    match mode:
        case "ppm":
            mz_diff = abs(mz2 - mz1) * 10^6 / mz2
        case "Da":
            mz_diff = abs(mz2 - mz1)
        case _:
            pass
    if mz_diff < tolerance:
        mz_score = 1 - mz_diff / tolerance
    else:
        mz_score = 0.0
    return mz_score


def get_spec_from_feature(feature: dict, 
                          spectra: dict):
    feature_id = feature["id"]
    return spectra.get(feature_id, None)
    

def get_ms2_spec_similarity(feature1: dict, 
                            feature2: dict, 
                            spectra: dict, 
                            mz_tol: float = 0.5, 
                            score_type: str = "None"):
    similarity_score = get_similarity_type(score_type, mz_tol)
    spec1 = get_spec_from_feature(feature1, spectra)
    spec2 = get_spec_from_feature(feature2, spectra)
    if spec1 is None or spec2 is None:
        return 0.0
    else:
        score = similarity_score.pair(spec1, spec2)
        return score['score'], score['matches']


def get_unannoted_feature(df: pd.DataFrame):
    return df[df['compound_name'] == '']


def get_annoted_feature(df1: pd.DataFrame, 
                        df2: pd.DataFrame):
    """
    Args:
        df1: query
        df2: reference
    """
    df = pd.merge(df1, df2, how='inner', on=['compound_name', 'adduct'])
    df = df.drop_duplicates(['id'])
    df = df[df['compound_name'] != '']
    return df 


def get_deep_lipid_name(annoted_feature: dict, 
                        unannoted_features: List[Dict], 
                        references: List[Dict], 
                        spectra: List[Spectrum],
                        rt_tol: float = 0.5, 
                        rt_mode: str = "absolute",
                        mz_tol: float = 1.0, 
                        mz_mode: str = "Da",
                        score_type: str = "None"):
    keys = ["id", "rt", "mz", "compound_name", "adduct", "cosine_score", "mol_formula", "smiles"]
    #annoted_feature = {key: annoted_feature[key] for key in keys}
    weight_mz = 0.25
    weight_rt = 0.25
    weight_spec = 0.5
    
    precursor_mz = annoted_feature["mz"]
    rt = annoted_feature["rt"]
    compound = annoted_feature["chain_shorthand"]
    

    unannoted_features = [filter_by_rt(feature, rt, rt_tol, rt_mode) for feature in unannoted_features]
    unannoted_features = [feature for feature in unannoted_features if feature is not None]
    
    candidates = [reference  for reference in references if reference["chain_shorthand"] == compound]
    
    if len(candidates) == 0:
        return compound
    
    candidate_score_list = []
    
    for candidate in candidates:
        G = nx.Graph()
        G.add_node(annoted_feature["id"], attributes = annoted_feature)
        
        delta_mass_list = candidate["delta_mass"]
        selected_feature_list = []
        for delta_mass in delta_mass_list:
            thoery_mz = precursor_mz + delta_mass
            selected_features = [filter_by_mz(feature, thoery_mz, mz_tol, mz_mode) for feature in unannoted_features]
            selected_features = [feature for feature in selected_features if feature is not None]
            
            selected_feature_list.extend(selected_features)
            
            if len(selected_features) > 0:
                for selected_feature in selected_features:
                    score_mz = get_precursor_mz_similarity(selected_feature, annoted_feature, mz_tol, mz_mode)
                    score_rt = get_rt_similarity(selected_feature, annoted_feature, rt_tol, rt_mode)
                    score_spec = get_ms2_spec_similarity(selected_feature, annoted_feature, spectra, mz_tol, score_type)
                    score = weight_mz * score_mz + weight_rt * score_rt + weight_spec * score_spec
                    selected_feature = {key: selected_feature[key] for key in keys}
                    selected_feature["score"] = score
                    G.add_node(selected_feature["id"], attributes = selected_feature)
                    G.add_edge(annoted_feature["id"], selected_feature["id"])                    
        candidate_score = 0
        for neighbor in G.neighbors(annoted_feature["id"]):   
            candidate_score += G.nodes[neighbor]["attributes"]["score"]
        #candidate_score /= len(delta_mass_list)
        candidate_score /= len(delta_mass_list)
        G.clear()
        #candidate_score_list.append(selected_feature_list)
        candidate_score_list.append(candidate_score)
    
    max_score = max(candidate_score_list)
    max_score_idx = candidate_score_list.index(max_score)
    lipid_name = candidates[max_score_idx]["lipid_name"]
    return lipid_name, max_score


def filter_by_rt(feature: dict, 
                 rt_theory: float, 
                 rt_tol: float = None, 
                 rt_mode: str = "absolute"):
    """filter by Retention Time
    """
    rt1 = feature["rt"]
    rt2 = rt_theory
    match rt_mode:
        case "relative":
            rt_diff = abs(rt1 - rt2) * 100 / rt2
        case "absolute":
            rt_diff = abs(rt1 - rt2)
        case _:
            pass
    if rt_diff < rt_tol:
        return feature
    else:
        return None
        

def filter_by_mz(feature: dict, 
                 mz_theory: float, 
                 mz_tol: float = None, 
                 mz_mode: str = "Da"):
    mz1 = feature["mz"]
    mz2 = mz_theory
    match mz_mode:
        case "ppm":
            mz_diff = abs(mz2 - mz1) * 10^6 / mz2
        case "Da":
            mz_diff = abs(mz2 - mz1)
        case _:
            pass
    if mz_diff < mz_tol:
        return feature
    else:
        return None
    

def algorithm_whole_pipeline(ms1_file,
                             ms2_file,
                             structure_file,
                             delta_mass_file,
                             out_file):
    #ms1_file = "data/processed_data/20240625-Hela/pos-O3/pos-Hela-O3_quant_full.csv"
    #ms2_file = "data/processed_data/20240625-Hela/pos-O3/pos-Hela-O3.mgf"
    #structure_file = "database/clean_data/reference_ms2_pos.tsv"
    #delta_mass_file = "database/clean_data/ozonolysis_delta_mass.json"
    #out_file = "results/20240702/pos-Hela-O3-ms1.csv"

    print("Running...")
    
    from .parse_utils import (
        load_structure_database,
        load_ozonolysis_prodcut_database,
        load_reference_database,
        load_ms1_peak_table,
        load_spec_file,
        
    )
    ## 0. Load reference name file & get ozID product delta mass
    df_structre = load_structure_database(structure_file)
    df_delta_mass = load_ozonolysis_prodcut_database(delta_mass_file)
    df_reference = load_reference_database(structure_file, delta_mass_file)

    ## 1. Load MS1 table
    df_ms1 = load_ms1_peak_table(ms1_file)

    ## 2. Get annoted table
    df_annoted  = get_annoted_feature(df_ms1, df_reference)
    df_unannoted = get_unannoted_feature(df_ms1)

    ## 3. Load MS/MS Spectra
    spectra = load_spec_file(ms2_file)
    
    ## 4.algorithm
    annoted_features = df_annoted.to_dict("records")
    unannoted_features = df_unannoted.to_dict("records")
    references = df_reference.to_dict("records")

    mix_list = [get_deep_lipid_name(annoted_feature, unannoted_features, references, spectra) for annoted_feature in annoted_features]
    lipid_name_list = [mix[0] for mix in mix_list]
    score_list = [mix[1] for mix in mix_list]

    df_final = df_annoted
    df_final = df_final.drop(["delta_mass"], axis=1)

    df_final['name'] = lipid_name_list
    df_final["score"] = score_list
    df_final.to_csv(out_file, index=None)
    
from typing import List, Tuple


import pandas as pd
import networkx as nx
from matchms import Spectrum


from .parse_utils import (
    load_structure_database,
    load_ozonolysis_prodcut_database,
    load_reference_database,
    load_ms1_peak_table,
    load_spec_file,
)

from .spectrum_utils import (
    get_rt_similarity,
    get_precursor_mz_similarity,
    get_ms2_spec_similarity,
    filter_by_rt,
    filter_by_mz
)


def get_unannoted_feature(df: pd.DataFrame) -> pd.DataFrame:
    """Get unannotated feature.

    Parameters:
        df:
            Query dataframe.

    Returns:
        `DataFrame` of unannotated feature
    """
    return df[df['compound_name'] == '']


def get_annoted_feature(
        df_query: pd.DataFrame,
        df_reference: pd.DataFrame
    ) -> pd.DataFrame:
    """Get unannotated feature.

    Parameters:
        df_query:
            Query dataframe.
        df_reference:
            Reference dataframe.

    Returns:
        df:
            `DataFrame` of annotated feature
    """
    df = pd.merge(df_query, df_reference, how='inner', on=['compound_name', 'adduct'])
    df = df.drop_duplicates(['id'])
    df = df[df['compound_name'] != '']

    return df


def get_deep_lipid_name(
        annoted_feature: dict,
        unannoted_features: List[dict],
        references: List[dict],
        spectra: List[Spectrum],
        rt_tol: float = 0.5,
        rt_mode: str = "absolute",
        mz_tol: float = 1.0,
        mz_mode: str = "Da",
        score_type: str = "None"
    ) -> Tuple[str, float]:
    """Lipid C=C identification.

    Parameters:
        annoted_feature:
            Annotated LC-MS feature.
        unannoted_features:
            Unannotated LC-MS features.
        references:
            Reference reference database.
        spectra:
            All of MS/MS spectra.
        rt_tol:
            Tolerance of retention time (RT).
        rt_mode:
            Type of retention time (RT) tolerance.
        mz_tol:
            Tolerance of m/z.
        mz_mode:
            Type of m/z tolerance.
        score_type:
            Type of MS/MS similarity score.

    Returns:
        lipid_name:
            Lipid with C=C information.
        max_score:
            Identification score.
    """
    keys = ["id", "rt", "mz", "compound_name", "adduct", "cosine_score", "mol_formula", "smiles"]
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
        candidate_score /= len(delta_mass_list)
        G.clear()
        candidate_score_list.append(candidate_score)

    max_score = max(candidate_score_list)
    max_score_idx = candidate_score_list.index(max_score)
    lipid_name = candidates[max_score_idx]["lipid_name"]

    return lipid_name, max_score


def algorithm_whole_pipeline(
        ms1_file: str,
        ms2_file: str,
        structure_file: str,
        delta_mass_file: str,
        out_file: str
    ):
    """Whole pipeline of SILICO-MS.

    Parameters:
        ms1_file:
            Filename of MS1 peak table.
        ms_file:
            Filename of MS2 spectra.
        structure_file:
            Filename of structure database.
        delta_mass_file:
            Filename of ozonolysis product database.
        out_file:
            Filename of output.

    Returns:
        None
    """

    print("Running...")

    ## 0. Load reference name file & get ozID product delta mass
    #df_structre = load_structure_database(structure_file)
    #df_delta_mass = load_ozonolysis_prodcut_database(delta_mass_file)
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

    mix_list = [get_deep_lipid_name(
                    annoted_feature=annoted_feature,
                    unannoted_features=unannoted_features,
                    references=references,
                    spectra=spectra
                ) for annoted_feature in annoted_features]
    lipid_name_list = [mix[0] for mix in mix_list]
    score_list = [mix[1] for mix in mix_list]

    df_final = df_annoted
    df_final = df_final.drop(["delta_mass"], axis=1)

    df_final['name'] = lipid_name_list
    df_final["score"] = score_list
    df_final.to_csv(out_file, index=None)

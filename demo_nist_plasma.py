"""Post-processing and main execution script for lipid annotation results.

This module provides utility functions for result filtering, lipid class counting,
YAML configuration handling, and the main pipeline execution. It processes raw
annotation outputs to generate clean, filtered lipid identifications and
statistical summaries of lipid class distributions.

Key Functions:
    dict_to_yaml: Save dictionary to YAML config file.
    yaml_to_dict: Load YAML config into dictionary.
    filter_results: Filter annotations by score and remove duplicates.
    count_lipid_class: Count and categorize lipids by lipid class.

Main Workflow:
    1. Load/ save configuration and parameters
    2. Run full lipid annotation pipeline
    3. Filter low-quality results
    4. Count lipid classes and export statistics
"""

import os
import yaml
import pandas as pd

from silico_ms.pipeline_runner import LipidAnnotatorRunner


def dict_to_yaml(
    file: str,
    data: dict
) -> None:
    """Write a dictionary to a YAML configuration file.

    Args:
        file: Path to output YAML file.
        data: Dictionary to be saved.
    """
    with open(file, "w", encoding="utf-8") as f:
        yaml.dump(data, f, allow_unicode=True, sort_keys=False)
    

def yaml_to_dict(file: str) -> dict:
    """Load a YAML configuration file into a Python dictionary.

    Args:
        file: Path to input YAML file.

    Returns:
        dict: Loaded configuration dictionary.
    """
    with open(file, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return config


def filter_results(
    data: pd.DataFrame,
    score_threshold: float = 0.1
) -> pd.DataFrame:
    """Filter lipid annotation results by quality and remove duplicates.

    Removes entries with total_score below threshold, filters out zero-intensity
    entries, and keeps only unique lipid structures per annotation.

    Args:
        data: Raw lipid annotation DataFrame.
        score_threshold: Minimum total score for valid annotations.

    Returns:
        pd.DataFrame: Filtered, high-quality lipid annotation results.
    """
    print("Filter by total score ...")
    
    var_cols = [
        "feature_id",	
        "rt",
        "mz_measured",
        "mz_calculated",
        "mz_error(ppm)",
        "mol_formula",
        "smiles",
        "adduct",
        "primary_annotated_name",
        "structure_annotated_name",
        "rt_score",
        "mz_score",
        "spec_score",
        "total_score",
        "fragments_ion_detection_level",
        "fragments_spectrum_detection_level",
    ]
    df = data.copy()
    df = df[df["total_score"] > score_threshold]
    print(df.shape)
    sample_cols = df.columns.difference(var_cols)
    df_filtered = df[df[sample_cols].fillna(0).sum(axis=1) > 0]
    print(df_filtered.shape)
    df_filtered = df_filtered.drop_duplicates(subset=["structure_annotated_name"])
        
    print("Filter by total score done!")
    
    return df_filtered


def count_lipid_class(
    data: pd.DataFrame
) -> pd.DataFrame:
    """Count identified lipids by their lipid class (PC, PE, TG, etc.).

    Maps lipid names to standard classes, counts occurrences,
    and outputs sorted results with a total row.

    Args:
        data: Filtered lipid annotation results.

    Returns:
        pd.DataFrame: Lipid class count statistics.
    """
    print("Count lipid class ...")
    
    mapping = {
        "PC": "PC",
        "PE": "PE",
        "PI": "PI",
        "PG": "PG",
        "PS": "PS",
        "DG": "DG",
        "TG": "TG"
    }
    specified_order = ["PC", "PE", "PI", "PG", "PS","DG", "TG", "Others"]
    specified_order = specified_order[::-1]

    df = data.copy()
    df["lipid_class"] = df["structure_annotated_name"].str.split(" ").str[0]
    df["lipid_class"] = df["lipid_class"].map(mapping).fillna("Others") # 其余的均记为"Others" 
    
    
    df_results = df["lipid_class"].value_counts()
    df_results = df_results.reset_index()
    df_results.columns = ["lipid_class", "count"]
    df_results["lipid_class"] = pd.Categorical(
                                    df_results["lipid_class"], 
                                    categories=specified_order, 
                                    ordered=True
                                )
    df_results = df_results.sort_values(by="lipid_class")
    total_count = df_results["count"].sum()
    df_results.loc[len(df_results)] = ["Total", total_count]
    
    print("Count lipid class done!")
    
    return df_results


if __name__ == "__main__":
    """Main execution entry for lipid annotation and post-processing pipeline."""
    
    # Configuration for NIST plasma dataset (positive + negative mode)
    NIST_PLASMA_CONFIG = {
        "pos":{
            "ozid_database_file": "database/clean_data/ozonolysis_delta_mass.json",
            "ozid": {
                "ms1_file": "example/data/nist_plasma/pos-nist_plasma-O3_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/nist_plasma/pos-nist_plasma-O3.mgf",
                "ms2_file_type": "mgf",
            },
            "cid": {
                "ms1_file": "example/data/nist_plasma/pos-nist_plasma-N2_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/nist_plasma/pos-nist_plasma-N2.mgf",
                "ms2_file_type": "mgf",
            }
        },
        "neg": {
            "ozid_database_file": "database/clean_data/ozonolysis_delta_mass.json",
            "ozid": {
                "ms1_file": "example/data/nist_plasma/neg-nist_plasma-O3_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/nist_plasma/neg-nist_plasma-O3.mgf",
                "ms2_file_type": "mgf",
            },
            "cid": {
                "ms1_file": "example/data/nist_plasma/neg-nist_plasma-N2_quant_full.csv",
                "ms1_file_type": "mzmine",
                "ms2_file": "example/data/nist_plasma/neg-nist_plasma-N2.mgf",
                "ms2_file_type": "mgf",
            }
        },
    }
    
    
    NIST_PLASMA_PARAMS = {
        "rt_tol": 0.5,
        "rt_tol_mode": "absolute",
        "mz_tol": 1.0,
        "mz_tol_mode": "Da",
        "ms2_spectrum_similarity_type": "ModifiedCosine",
        "score_threshold": 0.0,
        "sn_threshold": 2.0,
        "top_n": 10,
        "remove_false_positive": True,
        "base_out_dir": "example/hyperparams_opt/Hela",
        "rt_weight": 0.1,
        "precursor_mz_weight": 0.1,
        "spec_weight": 0.8,
        "out_dir": "example/results/20260109"
    }
    
    config_file = "example/data/nist_plasma/nist_plasma_config.yaml"
    params_file = "example/data/nist_plasma/nist_plasma_params_config.yaml"
    
    dict_to_yaml(
        file=config_file,
        data=NIST_PLASMA_CONFIG
    )
    dict_to_yaml(
        file=params_file,
        data=NIST_PLASMA_PARAMS
    )
    
    config = yaml_to_dict(file=config_file)
    params = yaml_to_dict(file=params_file)
    
    lipid_annotator_runner = LipidAnnotatorRunner(config=config, params=params)
    df_results_total = lipid_annotator_runner.annotate_lipids()
    df_filtered = filter_results(data=df_results_total, score_threshold=0.1)
    df_results = count_lipid_class(data=df_filtered)
    
    out_file = os.path.join(params["out_dir"], "lipid_class_count.csv")
    out_file = out_file.replace("\\", "/")
    print(out_file)
    df_results.to_csv(out_file, index=None)
    
"""Standalone lipid annotation script for mouse liver dataset.

This script runs end-to-end lipid annotation for positive and negative
ion mode LC-MS/MS data. It loads MS data, initializes annotation parameters,
runs lipid identification, performs false positive filtering,
and exports annotated results to CSV files.

Workflow:
    1. Load lipid fragment reference database
    2. Initialize feature processor and lipid annotator
    3. Load positive/negative ion mode MS1 and MS2 data
    4. Run lipid annotation for each ion mode
    5. Export annotation results and full feature tables
"""

from typing import Tuple, List
from matchms import Spectrum
import pandas as pd

from silico_ms.data_utils import DatabaseLoader, DataLoader
from silico_ms.spectrum_utils import (
    SpectrumFeatureProcesser, split_features, 
    ms2_spectra_to_dict
)
from silico_ms.algorithm import LipidAnnotator


def get_features(
    ms1_file: str,
    ms1_file_type: str,
    ms2_file: str,
    ms2_file_type: str
) -> Tuple[pd.DataFrame, pd.DataFrame, List[Spectrum]]:
    """Load and preprocess LC-MS/MS data into candidate and auxiliary features.

    Args:
        ms1_file: Path to MS1 peak table CSV file
        ms1_file_type: Format of MS1 file (e.g., mzmine)
        ms2_file: Path to MS2 spectrum file (e.g., mgf)
        ms2_file_type: Format of MS2 file (e.g., mgf)

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, List[Spectrum]]:
            Candidate features (annotated), auxiliary features (unannotated),
            and list of MS2 Spectrum objects
    """
    data_loader = DataLoader(
                    ms1_file=ms1_file, 
                    ms1_file_type=ms1_file_type, 
                    ms2_file=ms2_file, 
                    ms2_file_type=ms2_file_type
                )
    
    ms1_peak_table = data_loader.load_ms1_peak_table()
    ms2_spectra = data_loader.load_ms2_spectrum()
    df_candidate_features, df_auxiliary_features = split_features(
                                                ms1_peak_table=ms1_peak_table
                                            )
    return df_candidate_features, df_auxiliary_features, ms2_spectra


if __name__ == "__main__":
    ozid_database_file = "database/clean_data/ozonolysis_delta_mass.json"

    rt_tol = 0.5
    rt_tol_mode = "absolute"
    mz_tol = 1.0
    mz_tol_mode = "Da"
    rt_weight = 0.1
    precursor_mz_weight = 0.1
    spec_weight = 0.8
    ms2_spectrum_similarity_type = "None"
    score_threshold = 0.1
    sn_threshold = 2.0

    print("Load database")
    database_loader = DatabaseLoader(ozid_database_file=ozid_database_file)
    df_reference_database = database_loader.load_reference_database()
    print("Done!")
    
    feature_processer = SpectrumFeatureProcesser(
                            rt_tol=rt_tol,
                            rt_tol_mode=rt_tol_mode,
                            mz_tol=mz_tol,
                            mz_tol_mode=mz_tol_mode,
                            rt_weight=rt_weight,
                            precursor_mz_weight=precursor_mz_weight,
                            spec_weight=spec_weight,
                            ms2_spectrum_similarity_type=ms2_spectrum_similarity_type
                        )
    lipid_annotator = LipidAnnotator(
                        df_reference_database=df_reference_database, 
                        feature_processer=feature_processer, 
                        score_threshold=score_threshold,
                        top_n=3,
                        remove_false_positive=True,
                        sn_threshold=sn_threshold
                    )
    
    # pos mode
    pos_ms1_file = "example/data/mouse_liver1/pos_liver1_O3_quant_full.csv"
    pos_ms1_file_type = "mzmine"
    pos_ms2_file = "example/data/mouse_liver1/pos_liver1_O3.mgf"
    pos_ms2_file_type = "mgf"
    pos_ms1_out_file = "example/results/mouse_liver1/pos_liver1-ms1.csv"
    pos_features_all_out_file = "example/results/mouse_liver1/pos_liver1-ms1_all.csv"

    pos_ms1_file_blank = "example/data/mouse_liver1/pos_liver1_N2_quant_full.csv"
    pos_ms2_file_blank = "example/data/mouse_liver1/pos_liver1_N2.mgf"
    
    print("Load pos mode data")
    pos_df_candidate_features, pos_df_auxiliary_features, pos_ms2_spectra = get_features(
                                                        ms1_file=pos_ms1_file,
                                                        ms1_file_type=pos_ms1_file_type,
                                                        ms2_file=pos_ms2_file,
                                                        ms2_file_type=pos_ms2_file_type
                                                    )
    pos_ms2_spectra_dict = ms2_spectra_to_dict(ms2_spectra=pos_ms2_spectra)
    print("Done!")
    
    pos_data_loader_blank = DataLoader(
                                ms1_file=pos_ms1_file_blank, 
                                ms1_file_type=pos_ms1_file_type,
                                ms2_file=pos_ms2_file_blank, 
                                ms2_file_type=pos_ms2_file_type
                            )
    pos_df_features_blank = pos_data_loader_blank.load_ms1_peak_table()
    
    print("Process pos mode data")
    pos_df_results = lipid_annotator.annotate_features(
        df_candidate_features=pos_df_candidate_features, 
        df_auxiliary_features=pos_df_auxiliary_features,
        df_features_blank=pos_df_features_blank,
        ms2_spectra_dict=pos_ms2_spectra_dict,
    )
    pos_df_features_all = lipid_annotator.get_features_all()
    if pos_df_results is not None:
        pos_df_results.to_csv(pos_ms1_out_file, index=None)
        pos_df_features_all.to_csv(pos_features_all_out_file, index=None)
  
        print("Output Done!")

    print("Done!")
        
    # neg mode
    print("Process neg mode data")
    neg_ms1_file = "example/data/mouse_liver1/neg_liver1_O3_quant_full.csv"
    neg_ms1_file_type = "mzmine"
    neg_ms2_file = "example/data/mouse_liver1/neg_liver1_O3.mgf"
    neg_ms2_file_type = "mgf"

    neg_ms1_out_file = "example/results/mouse_liver1/neg_liver1-ms1.csv"
    neg_features_all_out_file = "example/results/mouse_liver1/neg_liver1-ms1_all.csv"
    
    neg_ms1_file_blank = "example/data/mouse_liver1/neg_liver1_N2_quant_full.csv"
    neg_ms2_file_blank = "example/data/mouse_liver1/neg_liver1_N2.mgf"
    
    print("Load neg mode data")
    neg_df_candidate_features, neg_df_auxiliary_features, neg_ms2_spectra = get_features(
                                                        ms1_file=neg_ms1_file,
                                                        ms1_file_type=neg_ms1_file_type,
                                                        ms2_file=neg_ms2_file,
                                                        ms2_file_type=neg_ms2_file_type
                                                    )
    neg_ms2_spectra_dict = ms2_spectra_to_dict(ms2_spectra=neg_ms2_spectra)
    print("Done!")
    
    neg_data_loader = DataLoader(
                        ms1_file=neg_ms1_file_blank, 
                        ms1_file_type=neg_ms1_file_type,
                        ms2_file=neg_ms2_file_blank, 
                        ms2_file_type=neg_ms2_file_type
                    )
    neg_df_features_blank = neg_data_loader.load_ms1_peak_table()
    
    neg_df_results = lipid_annotator.annotate_features(
        df_candidate_features=neg_df_candidate_features, 
        df_auxiliary_features=neg_df_auxiliary_features,
        df_features_blank=neg_df_features_blank,
        ms2_spectra_dict=neg_ms2_spectra_dict,
    )
    neg_df_features_all = lipid_annotator.get_features_all()
    
    if neg_df_results is not None:
        neg_df_results.to_csv(neg_ms1_out_file, index=None) 
        neg_df_features_all.to_csv(neg_features_all_out_file, index=None)
        print("Output Done!")
        
    print("Done!")
    
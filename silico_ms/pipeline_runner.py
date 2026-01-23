import os
from typing import Tuple


import pandas as pd
import yaml


from silico_ms.data_utils import DatabaseLoader, DataLoader
from silico_ms.spectrum_utils import (
    SpectrumFeatureProcesser, split_features, 
    ms2_spectra_to_dict
)
from silico_ms.algorithm import LipidAnnotator



class LipidAnnotatorRunner:
    def __init__(
        self,
        config: dict,
        params: dict
    ) -> None:
        self.config = config
        self.params = params
    
    def annotate_lipids(
        self,
    ) -> pd.DataFrame:
        """Annotate lipids with given hyperparameters."""
        out_dir = self.params["out_dir"]
        if not os.path.isdir(out_dir): 
            os.makedirs(out_dir)
        config_pos = self.config["pos"]
        config_neg = self.config["neg"]
        
        #print("Annotate lipids with hyperparameters:", out_dir)
        
        # pos
        print("Processing pos mode data ...")
        df_results_pos, df_results_all_pos = self._annotate_lipid_single(
                                    config=config_pos,
                                    params=self.params
                                )
        print("Done!")
        # neg
        print("Processing neg mode data ...")
        df_results_neg, df_results_all_neg = self._annotate_lipid_single(
                                    config=config_neg,
                                    params=self.params
                                )
        print("Done!")
        # total
        df_results_total = self._concat_df_results(
                                df_results_pos=df_results_pos, 
                                df_results_neg=df_results_neg
                            )
        df_results_all_total = self._concat_df_results(
                                df_results_pos=df_results_all_pos, 
                                df_results_neg=df_results_all_neg
                            )
        
        if df_results_pos is not None:
            out_file_pos = os.path.join(out_dir, "pos_feature_table.csv")
            out_file_all_pos = os.path.join(out_dir, "pos_feature_table_all.csv")
            self._write_dataframe(df=df_results_pos, out_file=out_file_pos)
            self._write_dataframe(df=df_results_all_pos, out_file=out_file_all_pos)
    
        if df_results_neg is not None:
            out_file_neg = os.path.join(out_dir, "neg_feature_table.csv")
            out_file_all_neg = os.path.join(out_dir, "neg_feature_table_all.csv")
            self._write_dataframe(df=df_results_neg, out_file=out_file_neg)
            self._write_dataframe(df=df_results_all_neg, out_file=out_file_all_neg)

        if df_results_total is not None:
            out_file_total = os.path.join(out_dir, "total_feature_table.csv")
            out_file_all_total = os.path.join(out_dir, "total_feature_table_all.csv")
            self._write_dataframe(df=df_results_total, out_file=out_file_total)
            self._write_dataframe(df=df_results_all_total, out_file=out_file_all_total)
        
        out_file_params = os.path.join(out_dir, "params.yaml")
        out_file_config = os.path.join(out_dir, "config.yaml")
        self._dict_to_yaml(data=self.config, file=out_file_config)
        self._dict_to_yaml(data=self.params, file=out_file_params)

        return df_results_total
    
    def _annotate_lipid_single(
        self,
        config: dict,
        params: dict
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        """
        ozid_database_file = config["ozid_database_file"]
        
        df_reference_database = self._load_database(ozid_database_file=ozid_database_file)
        lipid_annotator = self._set_params(
                                params=params,
                                df_reference_database=df_reference_database
                            )
        (df_candidate_features, 
         df_auxiliary_features, 
         df_features_cid, 
         ms2_spectra_dict) = self._load_data(config=config)


        df_results, df_results_all = self._annotate_features(
                                        lipid_annotator=lipid_annotator,
                                        df_candidate_features=df_candidate_features, 
                                        df_auxiliary_features=df_auxiliary_features,
                                        df_features_blank=df_features_cid,
                                        ms2_spectra_dict=ms2_spectra_dict
                                    )
        
              
        return df_results, df_results_all
        
    def _load_database(
        self,
        ozid_database_file: str
    ) -> pd.DataFrame:
        """
        """
        print("Loading database ...")
        database_loader = DatabaseLoader(ozid_database_file=ozid_database_file)
        df_reference_database = database_loader.load_reference_database()
        #print("Done!")
        
        return df_reference_database
    
    def _set_params(
        self,
        params: dict,
        df_reference_database: pd.DataFrame
    ) -> LipidAnnotator:
        """
        """
        rt_tol = params["rt_tol"]
        rt_tol_mode = params["rt_tol_mode"]
        mz_tol = params["mz_tol"]
        mz_tol_mode = params["mz_tol_mode"]
        
        rt_weight = params["rt_weight"]
        precursor_mz_weight = params["precursor_mz_weight"]
        spec_weight = params["spec_weight"]
        
        ms2_spectrum_similarity_type = params["ms2_spectrum_similarity_type"]
        score_threshold = params["score_threshold"]
        top_n = params["top_n"]
        remove_false_positive = params["remove_false_positive"]
        sn_threshold = params["sn_threshold"]

        
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
                            df_reference_database = df_reference_database, 
                            feature_processer = feature_processer, 
                            score_threshold = score_threshold,
                            top_n = top_n,
                            remove_false_positive = remove_false_positive,
                            sn_threshold = sn_threshold
                        )
        return lipid_annotator
    
    def _load_data(
        self,
        config: dict
    ) -> Tuple[pd.DataFrame]:
        """"""
        print("Loading data ...")
        ozid_ms1_file_ozid = config["ozid"]["ms1_file"]
        ozid_ms1_file_type = config["ozid"]["ms1_file_type"]
        ozid_ms2_file = config["ozid"]["ms2_file"]
        ozid_ms2_file_type = config["ozid"]["ms2_file_type"]
        cid_ms1_file = config["cid"]["ms1_file"]
        cid_ms1_file_type = config["cid"]["ms1_file_type"]
        cid_ms2_file = config["cid"]["ms2_file"]
        cid_ms2_file_type = config["cid"]["ms2_file_type"]
        
        data_loader_ozid = DataLoader(
                    ms1_file=ozid_ms1_file_ozid, 
                    ms1_file_type=ozid_ms1_file_type, 
                    ms2_file=ozid_ms2_file, 
                    ms2_file_type=ozid_ms2_file_type
                )
        data_loader_cid = DataLoader(
                                ms1_file=cid_ms1_file, 
                                ms1_file_type=cid_ms1_file_type,
                                ms2_file=cid_ms2_file, 
                                ms2_file_type=cid_ms2_file_type
                            )
        ms1_peak_table_ozid = data_loader_ozid.load_ms1_peak_table()
        ms2_spectra_ozid = data_loader_ozid.load_ms2_spectrum()
        ms2_spectra_dict = ms2_spectra_to_dict(ms2_spectra=ms2_spectra_ozid)
        df_candidate_features, df_auxiliary_features = split_features(
                                                ms1_peak_table=ms1_peak_table_ozid
                                            )
        df_features_cid = data_loader_cid.load_ms1_peak_table()
        #print("Done!")
        return df_candidate_features, df_auxiliary_features, df_features_cid, ms2_spectra_dict

    def _annotate_features(
        self,
        lipid_annotator: LipidAnnotator,
        df_candidate_features: pd.DataFrame,
        df_auxiliary_features: pd.DataFrame,
        df_features_blank: pd.DataFrame,
        ms2_spectra_dict: dict
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        """
        print("Annotate data ...")
        df_results = lipid_annotator.annotate_features(
            df_candidate_features=df_candidate_features, 
            df_auxiliary_features=df_auxiliary_features,
            df_features_blank=df_features_blank,
            ms2_spectra_dict=ms2_spectra_dict,
        )
        df_results_all = lipid_annotator.get_features_all()
        #print("Done!")
        
        return df_results, df_results_all
    
    def _concat_df_results(
        self,
        df_results_pos: pd.DataFrame,
        df_results_neg: pd.DataFrame
    ) -> pd.DataFrame:
        """
        """
        df_pos = df_results_pos.copy()
        df_neg = df_results_neg.copy()
        df_pos["feature_id"] = "P" + df_pos["feature_id"].astype(str)
        df_neg["feature_id"] = "N" + df_neg["feature_id"].astype(str)
        
        df_results_total = pd.concat(
                        [df_pos, df_neg], 
                        axis=0, 
                        ignore_index=True
                    )
        return df_results_total

    def _write_dataframe(
        self,
        df: pd.DataFrame,
        out_file: str
    ) -> None:
        """
        """
        if df is not None:
            df.to_csv(out_file, index=None)

    def _dict_to_yaml(
        self,
        data: dict,
        file: str
    ) -> None:
        """Write a dictionary to a YAML file."""
        with open(file, "w", encoding="utf-8") as f:
            yaml.dump(data, f, allow_unicode=True, sort_keys=False)

from typing import List
import pandas as pd
from tqdm import tqdm
tqdm.pandas(desc="Lipid C=C annotation") 


from silico_ms.spectrum_utils import SpectrumFeatureProcesser


class LipidAnnotator:
    def __init__(
        self,
        df_reference_database: pd.DataFrame,
        feature_processer: SpectrumFeatureProcesser,
        score_threshold: float = 0.0,
        top_n: int = 3,
        remove_false_positive: bool = False,
        sn_threshold: float = 1.0,
        epsilon: float = 1e-3
    ):
        self.reference_database = df_reference_database
        self.feature_processer = feature_processer
        self.score_threshold = score_threshold
        self.top_n = top_n
        self.remove_false_positive = remove_false_positive
        self.sn_threshold = sn_threshold
        self.epsilon = epsilon

    def annotate_features(
        self,
        df_candidate_features: pd.DataFrame,
        df_auxiliary_features: pd.DataFrame,
        df_features_blank: pd.DataFrame,
        ms2_spectra_dict: dict,
    ) -> pd.DataFrame:
        """"""
        self.features_all = pd.DataFrame()
        columns_to_drop = [
            "fragment_mz_delta",
            "fragment_fatty_acid_pos",
            "fragment_types"
        ]
    
        df_result_list = df_candidate_features.progress_apply(
            self._annotate_feature_single, 
            axis=1, 
            df_auxiliary_features=df_auxiliary_features,
            df_features_blank=df_features_blank,
            ms2_spectra_dict=ms2_spectra_dict,
        ).tolist()

        df_results = pd.concat(
            df_result_list,
            axis=0,
            ignore_index=True
        )
        df_results = df_results.drop(columns=columns_to_drop)
        
        return df_results

    def _annotate_feature_single(
        self,
        candidate_feature: pd.Series,
        df_auxiliary_features: pd.DataFrame,
        df_features_blank: pd.DataFrame,
        ms2_spectra_dict: dict
    ) -> pd.DataFrame:
        """"""
        df_reference_lipids = self.get_refence_lipids(
                                candidate_feature=candidate_feature
                            )
        
        if df_reference_lipids.empty:
            return pd.DataFrame()
        
        df_fragment_features = self.generate_df_fragment_features(
                                candidate_feature=candidate_feature,
                                df_auxiliary_features=df_auxiliary_features,
                                df_reference_lipids=df_reference_lipids,
                                ms2_spectra_dict=ms2_spectra_dict
                            )
        if self.remove_false_positive:
            if df_fragment_features.empty:
                return pd.DataFrame()
        
            df_fragment_features = self.filter_false_positive_fragment_features(
                                    candidate_feature=candidate_feature, 
                                    df_fragment_features=df_fragment_features,
                                    df_features_blank=df_features_blank
                                )
        
        if df_fragment_features.empty:
            return pd.DataFrame()
        
        df_result_features = self.generate_df_result_features(
                                    candidate_feature=candidate_feature,
                                    df_fragment_features=df_fragment_features
                                )
        df_result_features_filter = self.get_top_n_candidates(
                                df_result_features=df_result_features,
                            )
        
        
        if df_result_features_filter.empty:
            return pd.DataFrame()
        
        df_fragment_features_filter = self.get_top_n_fragments(
                                            df_result_features=df_result_features_filter,
                                            df_fragment_features=df_fragment_features
                                        )
        df_result = self.quantification(
                                df_result_features=df_result_features_filter,
                                df_fragment_features=df_fragment_features_filter
                            )
        
        self.features_all = pd.concat(
            objs=[
                self.features_all, 
                df_result, 
                df_fragment_features_filter
            ], 
            axis=0, 
            ignore_index=True
        )
        
        return df_result

    def get_refence_lipids(
        self,
        candidate_feature: dict
    ) -> pd.DataFrame:
        """"""
        primary_annotated_name = candidate_feature.get("primary_annotated_name")
        df = self.reference_database.copy()
        df_reference_lipids = df[df["primary_annotated_name"] == primary_annotated_name]
        
        return df_reference_lipids

    def generate_df_fragment_features(
        self,
        df_reference_lipids: pd.DataFrame,
        candidate_feature: pd.Series,
        df_auxiliary_features: pd.DataFrame,
        ms2_spectra_dict: dict
    ) -> pd.DataFrame:
        df_fragment_features_list = df_reference_lipids.apply(
            self.get_fragment_features,
            axis=1,
            candidate_feature=candidate_feature,
            df_auxiliary_features=df_auxiliary_features,
            ms2_spectra_dict=ms2_spectra_dict
        ).tolist()
        df_fragment_features = pd.concat(
                                df_fragment_features_list, 
                                axis=0,
                                ignore_index=True
                            )
        if df_fragment_features.empty:
            return pd.DataFrame()

        df_fragment_features = self.filter_fragment_features(
                                    df_fragment_features=df_fragment_features
                                )
        return df_fragment_features
    
    def get_fragment_features(
        self,
        reference_lipid: pd.Series,
        candidate_feature: pd.Series,
        df_auxiliary_features: pd.DataFrame,
        ms2_spectra_dict: dict
    ) -> pd.DataFrame:
        """"""
        df_features = self.feature_processer.filter_by_rt(
                                    df_features=df_auxiliary_features,
                                    rt_theory=candidate_feature.get("rt")
                                )
        
        df_fragment_features_list = []
        fragment_mz_delta = reference_lipid.get("fragment_mz_delta")
        fragment_fatty_acid_pos = reference_lipid.get("fragment_fatty_acid_pos")
        fragment_types = reference_lipid.get("fragment_types")

        for idx, mz_delta in enumerate(fragment_mz_delta):
            df_tmp_features = self.feature_processer.filter_by_precursor_mz(
                    df_features=df_features,
                    precursor_mz_theory=candidate_feature.get("precursor_mz") + mz_delta
                )
            if df_tmp_features.empty:
                continue
                
            fragment_scores = df_tmp_features.apply(
                self.feature_processer.feature_similarity,
                axis=1,
                feature2=candidate_feature,
                mz_delta=mz_delta,
                ms2_spectra_dict=ms2_spectra_dict
            )
            df_tmp_features["fragment_mz_delta"] = mz_delta
            df_tmp_features["primary_annotated_name"] = reference_lipid.get("primary_annotated_name")
            df_tmp_features["structure_annotated_name"] = reference_lipid.get("structure_annotated_name")
            df_tmp_features["fragment_fatty_acid_pos"] = fragment_fatty_acid_pos[idx]
            df_tmp_features["fragment_types"] = fragment_types[idx]
            df_tmp_features["score"] = fragment_scores
            
            df_fragment_features_list.append(df_tmp_features)
        
        if not df_fragment_features_list:
            return pd.DataFrame()

        df_fragment_features = pd.concat(
            df_fragment_features_list, 
            axis=0, 
            ignore_index=True
        )

        return df_fragment_features

    def filter_fragment_features(
        self,
        df_fragment_features: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        """
        idx = df_fragment_features.groupby(
           by=[
                "structure_annotated_name", 
                "fragment_fatty_acid_pos", 
                "fragment_types"
            ]
        )["score"].idxmax()
        df_fragment_features_filtered = df_fragment_features.loc[idx]
        
        return df_fragment_features_filtered
    
    
    def filter_false_positive_fragment_features(
        self,
        candidate_feature: pd.Series,
        df_fragment_features: pd.DataFrame,
        df_features_blank: pd.DataFrame
    ) -> pd.DataFrame:
        """"""
        candidate_feature_blank = self.feature_processer.search_similar_feature(
                                    feature=candidate_feature, 
                                    df_features_blank=df_features_blank
                                ).squeeze()
        if candidate_feature_blank.empty:
            return df_fragment_features.copy()
    
        df_fragment_features_blank_list = df_fragment_features.apply(
                                            self.feature_processer.search_similar_feature, 
                                            axis=1,
                                            df_features_blank=df_features_blank
                                        ).values
        
        if len(df_fragment_features_blank_list) == 0:
            return df_fragment_features.copy()
        
        df_fragment_features_blank = pd.concat(
                                df_fragment_features_blank_list, 
                                axis=0,
                                ignore_index=True
                            )
        if df_fragment_features_blank.empty:
            return df_fragment_features.copy()
        
        df_fragment_features_norm = self.normalize_features(
                                        df_fragment_features=df_fragment_features,
                                        candidate_feature=candidate_feature
                                    )

        df_fragment_features_blank_norm = self.normalize_features(
                                            df_fragment_features=df_fragment_features_blank,
                                            candidate_feature=candidate_feature_blank
                                        )
        df_fragment_features_filter = self.remove_false_positive_features(
                                            df_features=df_fragment_features_norm, 
                                            df_features_blank=df_fragment_features_blank_norm
                                        )

        if df_fragment_features_filter.empty:
            return pd.DataFrame()
        
        return df_fragment_features_filter
        
    def normalize_features(
        self,
        df_fragment_features: pd.DataFrame,
        candidate_feature: pd.Series
    ) -> pd.DataFrame:
        """"""
        df = df_fragment_features.copy()
        var_columns = [
            "feature_id",
            "rt",
            "precursor_mz",
            "adduct",
            "fragment_mz_delta",
            "primary_annotated_name",
            "structure_annotated_name",
            "fragment_fatty_acid_pos",
            "fragment_types",
            "score",
        ]
        sample_columns = df.columns.difference(var_columns).tolist()
        denom = candidate_feature[sample_columns].mask(lambda x: x == 0.0, self.epsilon)
        df[sample_columns] = df[sample_columns].div(
                                denom, 
                                axis=1
                            ) * 100
        
        return df
     
    def remove_false_positive_features(
        self,
        df_features: pd.DataFrame,
        df_features_blank: pd.DataFrame
    ) -> pd.DataFrame:
        """"""
        df1 = df_features.copy()
        df2 = df_features_blank.copy()
        
        key_cols = [
            "structure_annotated_name",
            "fragment_fatty_acid_pos",
            "fragment_types",
        ]
        var_columns = [
            "feature_id",
            "rt",
            "precursor_mz",
            "adduct",
            "fragment_mz_delta",
            "primary_annotated_name",
            "structure_annotated_name",
            "fragment_fatty_acid_pos",
            "fragment_types",
            "score",
        ]
        sample_columns = df1.columns.difference(var_columns).tolist()
        
        df1["key"] = df1[key_cols].astype(str).agg("||".join, axis=1)
        df2["key"] = df2[key_cols].astype(str).agg("||".join, axis=1)
        
        df_tmp = pd.merge(
            left=df1,
            right=df2[["key"] + sample_columns],
            on="key",
            how="left",
            suffixes=("", "_blank"),
            indicator=True
        )
        df_mask = pd.DataFrame(
            index=df_tmp.index, 
            columns=sample_columns, 
            dtype=float
        )
        
        for col in sample_columns:
            unmatched = (df_tmp["_merge"] == "left_only")
            denominator = df_tmp[f"{col}_blank"]
            denominator = denominator.mask(denominator.isna(), self.epsilon)
            denominator = denominator.mask(lambda x: x == 0.0, self.epsilon)
            ratio = df_tmp[col] / denominator
            matched_flag = (ratio >= self.sn_threshold)
            df_mask[col] = matched_flag.mask(unmatched, True)
        
        row_mask = df_mask.any(axis=1)
        df_result = df_tmp[row_mask]
        delete_columns = [
            c 
            for c in df_result.columns
            if str(c).endswith("_blank")
        ] + ["key", "_merge"]
        df_result = df_result.drop(columns=delete_columns)
        
        return df_result
        
    def generate_df_result_features(
        self,
        candidate_feature: pd.Series,
        df_fragment_features: pd.DataFrame
    ) -> pd.DataFrame:
        df_candidate_feature = candidate_feature.to_frame().T
        df_result_feature = pd.merge(
            left=df_candidate_feature, 
            right=df_fragment_features[[
                        "primary_annotated_name", 
                        "structure_annotated_name"
                    ]].drop_duplicates(),
            on=["primary_annotated_name"],
            how="inner"
        )
        
        result_scores = self.score_result_feature(df_fragment_features=df_fragment_features)
        df_result_feature["fragment_mz_delta"] = 0
        df_result_feature["fragment_fatty_acid_pos"] = 0
        df_result_feature["fragment_types"] = "substrate"
        df_result_feature["score"] = result_scores
        
        return df_result_feature
    
    def score_result_feature(
        self,
        df_fragment_features: pd.DataFrame,
    ) -> List[float]:
        result_scores = df_fragment_features.groupby(
            by=["structure_annotated_name"]
        )["score"].mean().tolist()
        
        return result_scores
        
    def get_top_n_candidates(
        self,
        df_result_features: pd.DataFrame,
    ) -> pd.DataFrame:
        """"""
        df_tmp = df_result_features.nlargest(n=self.top_n, columns="score")
        df_tmp = df_tmp[df_tmp["score"] > self.score_threshold]
        
        if df_tmp.empty:
            df_result_filter = df_result_features.nlargest(n=1, columns="score")
        else:
            df_result_filter = df_tmp
        
        return df_result_filter

    def get_top_n_fragments(
      self,
      df_result_features: pd.DataFrame,
      df_fragment_features: pd.DataFrame  
    ) -> pd.DataFrame:
        """"""
        df = df_fragment_features.copy()
        name = df_result_features["structure_annotated_name"].tolist()
        df_fragment_features_filter = df[df["structure_annotated_name"].isin(name)]
        
        return df_fragment_features_filter
    
    def quantification(
        self,
        df_result_features: pd.DataFrame,
        df_fragment_features: pd.DataFrame
    ) -> pd.DataFrame:
        df_result = df_result_features.copy()
        var_columns = [
            "feature_id",
            "rt",
            "precursor_mz",
            "adduct",
            "fragment_mz_delta",
            "primary_annotated_name",
            "structure_annotated_name",
            "fragment_fatty_acid_pos",
            "fragment_types",
            "score",
        ]
        sample_columns = df_fragment_features.columns.difference(var_columns).tolist()
        df_intnsities_sum = df_fragment_features.groupby(
            by=["structure_annotated_name"]
        )[sample_columns].sum().reset_index()
        
        df_intnsities_sum[sample_columns] = df_intnsities_sum[sample_columns].astype(float)
        totals = df_intnsities_sum[sample_columns].sum()
        df_intnsities_norm = df_intnsities_sum[sample_columns].div(totals)
        
        norm_factors = df_intnsities_norm[sample_columns].values
        df_result[sample_columns] *= norm_factors
        
        return df_result
        
    def get_features_all(self) -> pd.DataFrame:
        return self.features_all
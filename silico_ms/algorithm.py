"""Lipid annotation and scoring module for LC-MS/MS lipidomics data analysis.

This module implements the core lipid annotation logic using reference databases
and feature similarity matching. It identifies lipid C=C double bond positions
by matching candidate features with theoretical fragment patterns.

Key functionalities:
1. Annotate candidate lipid features using a reference lipid database
2. Filter and score fragment features by RT, m/z, and MS2 spectrum similarity
3. Remove false positive annotations using blank subtraction
4. Calculate confidence scores and detection levels for lipid identifications
5. Generate quantified and formatted annotation results

Classes:
    LipidAnnotator: Main class for automated lipid annotation and scoring.
"""

from typing import List
import pandas as pd
from tqdm import tqdm
tqdm.pandas(desc="Lipid C=C annotation") 

from silico_ms.spectrum_utils import SpectrumFeatureProcesser


class LipidAnnotator:
    """Main class for lipid feature annotation and scoring.

    Manages the annotation pipeline: matching candidate features to reference
    lipid fragments, filtering by similarity thresholds, removing false positives,
    and computing final annotation scores.

    Attributes:
        reference_database: DataFrame of reference lipid fragment patterns.
        feature_processer: Instance for feature filtering and similarity calculation.
        score_threshold: Minimum total score for valid annotations.
        top_n: Number of top-scoring annotations to retain.
        remove_false_positive: Enable blank subtraction for noise filtering.
        sn_threshold: S/N threshold for false positive filtering.
        epsilon: Small value to avoid division by zero.
        var_columns: List of metadata/annotation columns.
    """
    
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
        """Initialize LipidAnnotator with parameters and reference data.

        Args:
            df_reference_database: Reference lipid fragment database.
            feature_processer: Pre-configured feature processor.
            score_threshold: Minimum score threshold for annotations.
            top_n: Keep top N highest-scoring lipid structures.
            remove_false_positive: Enable blank-based noise removal.
            sn_threshold: Signal-to-noise threshold for filtering.
            epsilon: Small constant for numerical stability.
        """
        self.reference_database = df_reference_database
        self.feature_processer = feature_processer
        self.score_threshold = score_threshold
        self.top_n = top_n
        self.remove_false_positive = remove_false_positive
        self.sn_threshold = sn_threshold
        self.epsilon = epsilon
        self.var_columns = self.set_var_columns()
        
    def set_var_columns(self) -> List[str]:
        """Define fixed metadata and annotation column names.

        Returns:
            List[str]: List of core feature annotation columns.
        """
        var_columns = [
            "feature_id",
            "rt",
            "precursor_mz",
            "mz_calculated",
            "mol_formula",
            "smiles",
            "adduct",
            "fragment_mz_delta",
            "primary_annotated_name",
            "structure_annotated_name",
            "fragment_fatty_acid_pos",
            "fragment_types",
            "rt_score",
            "mz_score",
            "spec_score",
            "total_score",
            "fragments_ion_detection_level",
            "fragments_spectrum_detection_level",
        ]
        return var_columns

    def annotate_features(
        self,
        df_candidate_features: pd.DataFrame,
        df_auxiliary_features: pd.DataFrame,
        df_features_blank: pd.DataFrame,
        ms2_spectra_dict: dict,
    ) -> pd.DataFrame:
        """Main pipeline to annotate all candidate lipid features.

        Processes each candidate feature by matching against reference fragments,
        scoring, filtering, and aggregating results into a final DataFrame.

        Args:
            df_candidate_features: Annotated candidate features.
            df_auxiliary_features: Unannotated fragment candidates.
            df_features_blank: Blank/control features for noise removal.
            ms2_spectra_dict: MS2 spectra lookup by feature_id.

        Returns:
            pd.DataFrame: Final annotated lipid results.
        """
        self.features_all = pd.DataFrame()
        df_result_list = df_candidate_features.progress_apply(
            self._annotate_feature_single, 
            axis=1, 
            df_auxiliary_features=df_auxiliary_features,
            df_features_blank=df_features_blank,
            ms2_spectra_dict=ms2_spectra_dict,
        ).tolist()
        df_result = self.process_df_result_list(df_list=df_result_list)
        
        self.features_all = self.process_df_features_all(df=self.features_all)
        
        return df_result

    def _annotate_feature_single(
        self,
        candidate_feature: pd.Series,
        df_auxiliary_features: pd.DataFrame,
        df_features_blank: pd.DataFrame,
        ms2_spectra_dict: dict
    ) -> pd.DataFrame:
        """Annotate one single candidate feature.

        Internal workflow:
        1. Retrieve reference lipids
        2. Find matching fragment features
        3. Optional false positive filtering
        4. Score and select top candidates
        5. Quantify and assemble output

        Args:
            candidate_feature: Single lipid feature to annotate.
            df_auxiliary_features: Pool of potential fragment features.
            df_features_blank: Blank features for noise removal.
            ms2_spectra_dict: MS2 spectra lookup dict.

        Returns:
            pd.DataFrame: Annotation result for the candidate.
        """
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
        
        df_result = df_result.dropna(axis=1, how="all")
        df_fragment_features_filter = df_fragment_features_filter.dropna(axis=1, how="all")
        
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
        """Retrieve reference lipid entries matching the candidate annotation.

        Args:
            candidate_feature: Feature with primary_annotated_name.

        Returns:
            pd.DataFrame: Matching reference lipid fragments.
        """
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
        """Generate fragment features by matching reference patterns to observed features.

        Args:
            df_reference_lipids: Reference fragment patterns.
            candidate_feature: Parent lipid feature.
            df_auxiliary_features: Pool of potential fragments.
            ms2_spectra_dict: MS2 spectra lookup.

        Returns:
            pd.DataFrame: Scored fragment features.
        """
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

        return df_fragment_features
    
    def get_fragment_features(
        self,
        reference_lipid: pd.Series,
        candidate_feature: pd.Series,
        df_auxiliary_features: pd.DataFrame,
        ms2_spectra_dict: dict
    ) -> pd.DataFrame:
        """Find and score fragment features for one reference lipid structure.

        Args:
            reference_lipid: Single reference lipid entry.
            candidate_feature: Parent feature.
            df_auxiliary_features: Pool of potential fragments.
            ms2_spectra_dict: MS2 spectra lookup.

        Returns:
            pd.DataFrame: Scored fragment features for this structure.
        """
        df_features = self.feature_processer.filter_by_rt(
                                    df_features=df_auxiliary_features,
                                    rt_theory=candidate_feature.get("rt")
                                )
        
        df_fragment_features_list: List[pd.DataFrame] = []
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
            
            df_fragment_scores = df_tmp_features.apply(
                self.feature_processer.feature_similarity,
                axis=1,
                feature2=candidate_feature,
                mz_delta=mz_delta,
                ms2_spectra_dict=ms2_spectra_dict
            )
            if df_fragment_scores.empty:
                continue
            
            df_tmp_features["fragment_mz_delta"] = mz_delta
            df_tmp_features["primary_annotated_name"] = reference_lipid.get("primary_annotated_name")
            df_tmp_features["structure_annotated_name"] = reference_lipid.get("structure_annotated_name")
            df_tmp_features["fragment_fatty_acid_pos"] = fragment_fatty_acid_pos[idx]
            df_tmp_features["fragment_types"] = fragment_types[idx]
            df_tmp_features = df_tmp_features.reset_index(drop=True)
            df_fragment_scores = df_fragment_scores.reset_index(drop=True)
            df_tmp_features = pd.concat([df_tmp_features, df_fragment_scores], axis=1)
            # add smiles
            df_tmp_features["smiles"] = reference_lipid.get("smiles")
            
            df_fragment_features_list.append(df_tmp_features)
        
        if not df_fragment_features_list:
            return pd.DataFrame()

        df_fragment_features = pd.concat(
            df_fragment_features_list, 
            axis=0, 
            ignore_index=True
        )
        df_fragment_features = self.filter_fragment_features(df_fragment_features)
        
        max_fragments = len(fragment_mz_delta)
        fragments_ion_detection_level = self._get_fragment_ion_level(
                                df_fragment_features=df_fragment_features, 
                                max_fragments=max_fragments
                            )
        ms2_spectrum_level = self._get_fragment_ms_spec_level(
                                df_fragment_features=df_fragment_features, 
                                max_fragments=max_fragments, 
                                ms2_spectra_dict=ms2_spectra_dict
                            )
        df_fragment_features["fragments_ion_detection_level"] = fragments_ion_detection_level
        df_fragment_features["fragments_spectrum_detection_level"] = ms2_spectrum_level
        
        return df_fragment_features

    def filter_fragment_features(
        self,
        df_fragment_features: pd.DataFrame,
    ) -> pd.DataFrame:
        """Keep only the highest-scoring fragment per structure and position.

        Args:
            df_fragment_features: Scored fragment features.

        Returns:
            pd.DataFrame: Filtered best fragments.
        """
        idx = df_fragment_features.groupby(
           by=[
                "structure_annotated_name", 
                "fragment_fatty_acid_pos", 
                "fragment_types"
            ]
        )["total_score"].idxmax()
        df_fragment_features_filtered = df_fragment_features.loc[idx]
        
        return df_fragment_features_filtered
    
    def _get_fragment_ion_level(
        self, 
        df_fragment_features: pd.DataFrame, 
        max_fragments: int
    ) -> str:
        """Classify fragment ion detection level: None / Some / All.

        Args:
            df_fragment_features: Matched fragment features.
            max_fragments: Expected number of fragments.

        Returns:
            str: Detection level label.
        """
        num_fragments = len(df_fragment_features)
        if num_fragments == 0:
            fragments_ion_detection_level = "None"
        elif num_fragments < max_fragments:
            fragments_ion_detection_level = "Some"
        elif num_fragments == max_fragments:
             fragments_ion_detection_level = "All"
        else:
            raise ValueError(f"Max fragment is {max_fragments}, \
                             but get {num_fragments}")

        return fragments_ion_detection_level
    
    def _get_fragment_ms_spec_level(
        self, 
        df_fragment_features: pd.DataFrame, 
        max_fragments: int, 
        ms2_spectra_dict: dict
    ) -> str:
        """Classify MS2 spectrum detection level for fragments.

        Args:
            df_fragment_features: Matched fragment features.
            max_fragments: Expected fragments.
            ms2_spectra_dict: MS2 spectra lookup.

        Returns:
            str: Spectrum detection level.
        """
        df = df_fragment_features.copy()
        df_filtered = df[df["fragment_mz_delta"] < 0]
        
        num_fragments_spec = df_filtered["feature_id"].isin(ms2_spectra_dict.keys()).sum()
        if num_fragments_spec == 0:
            fragments_spec_level = "None"
        elif num_fragments_spec < max_fragments:
            fragments_spec_level = "Some"
        elif num_fragments_spec == max_fragments:
             fragments_spec_level = "All"
        else:
            raise ValueError(f"Max fragment is {max_fragments}, \
                             but get {num_fragments_spec}")
            
        return fragments_spec_level
    
    def filter_false_positive_fragment_features(
        self,
        candidate_feature: pd.Series,
        df_fragment_features: pd.DataFrame,
        df_features_blank: pd.DataFrame
    ) -> pd.DataFrame:
        """Remove false positive fragments using blank/control sample subtraction.

        Args:
            candidate_feature: Parent feature.
            df_fragment_features: Matched fragments.
            df_features_blank: Blank sample features.

        Returns:
            pd.DataFrame: Filtered fragments without false positives.
        """
    
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
        """Normalize fragment intensities relative to parent feature.

        Args:
            df_fragment_features: Fragment features.
            candidate_feature: Parent feature for normalization.

        Returns:
            pd.DataFrame: Normalized fragment features.
        """
        df = df_fragment_features.copy()
        sample_columns = df.columns.difference(self.var_columns).tolist()
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
        """Remove fragments present in blank sample above threshold.

        Args:
            df_features: Sample fragment features.
            df_features_blank: Blank fragment features.

        Returns:
            pd.DataFrame: Fragments passing blank subtraction.
        """
        df1 = df_features.copy()
        df2 = df_features_blank.copy()
        
        key_cols = [
            "structure_annotated_name",
            "fragment_fatty_acid_pos",
            "fragment_types",
        ]
        sample_columns = df1.columns.difference(self.var_columns).tolist()
        
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
        """Build final candidate annotation with averaged fragment scores.

        Args:
            candidate_feature: Parent feature.
            df_fragment_features: Scored fragments.

        Returns:
            pd.DataFrame: Candidate feature with annotation scores.
        """
        df_candidate_feature = candidate_feature.to_frame().T
        df_result_feature = pd.merge(
            left=df_candidate_feature, 
            right=df_fragment_features[[
                        "primary_annotated_name", 
                        "structure_annotated_name",
                        # add smiles
                        "smiles",
                        "fragments_ion_detection_level",
                        "fragments_spectrum_detection_level",
                    ]].drop_duplicates(),
            on=["primary_annotated_name"],
            how="inner"
        )
        
        df_result_scores = self.score_result_feature(df_fragment_features=df_fragment_features)
        df_result_feature["fragment_mz_delta"] = 0
        df_result_feature["fragment_fatty_acid_pos"] = 0
        df_result_feature["fragment_types"] = "substrate"
        df_result_feature = df_result_feature.reset_index(drop=True)
        df_result_scores = df_result_scores.reset_index(drop=True)
        df_result_feature = pd.concat([df_result_feature, df_result_scores], axis=1)
        
        return df_result_feature
    
    def score_result_feature(
        self,
        df_fragment_features: pd.DataFrame,
    ) -> pd.DataFrame:
        """Compute mean scores across all fragments for each structure.

        Args:
            df_fragment_features: Scored fragment features.

        Returns:
            pd.DataFrame: Mean scores per structure.
        """
        df_result_scores = df_fragment_features.groupby(
            by=["structure_annotated_name"]
        )[["rt_score","mz_score", "spec_score","total_score"]].mean()
        
        return df_result_scores
        
    def get_top_n_candidates(
        self,
        df_result_features: pd.DataFrame,
    ) -> pd.DataFrame:
        """Select top N highest-scoring annotations above threshold.

        Args:
            df_result_features: Scored candidate annotations.

        Returns:
            pd.DataFrame: Top N filtered results.
        """
        df_tmp = df_result_features.nlargest(n=self.top_n, columns="total_score")
        df_tmp = df_tmp[df_tmp["total_score"] > self.score_threshold]
        
        if df_tmp.empty:
            df_result_filter = df_result_features.nlargest(n=1, columns="total_score")
        else:
            df_result_filter = df_tmp
        
        return df_result_filter

    def get_top_n_fragments(
      self,
      df_result_features: pd.DataFrame,
      df_fragment_features: pd.DataFrame  
    ) -> pd.DataFrame:
        """Retain fragments corresponding to top N candidate structures.

        Args:
            df_result_features: Top N candidates.
            df_fragment_features: All scored fragments.

        Returns:
            pd.DataFrame: Fragments for top candidates.
        """
        df = df_fragment_features.copy()
        name = df_result_features["structure_annotated_name"].tolist()
        df_fragment_features_filter = df[df["structure_annotated_name"].isin(name)]
        
        return df_fragment_features_filter
    
    def quantification(
        self,
        df_result_features: pd.DataFrame,
        df_fragment_features: pd.DataFrame
    ) -> pd.DataFrame:
        """Quantify sample intensities using fragment sum normalization.

        Args:
            df_result_features: Top candidate annotations.
            df_fragment_features: Corresponding fragments.

        Returns:
            pd.DataFrame: Quantified annotation results.
        """
        df_result = df_result_features.copy()
        sample_columns = df_fragment_features.columns.difference(self.var_columns).tolist()
        df_intnsities_sum = df_fragment_features.groupby(
            by=["structure_annotated_name"]
        )[sample_columns].sum().reset_index()
        
        df_intnsities_sum[sample_columns] = df_intnsities_sum[sample_columns].astype(float)
        totals = df_intnsities_sum[sample_columns].sum()
        df_intnsities_norm = df_intnsities_sum[sample_columns].div(totals)
        
        norm_factors = df_intnsities_norm[sample_columns].values
        df_result[sample_columns] *= norm_factors
        
        return df_result
    
    def process_df_result_list(
        self,
        df_list: list
    ) -> pd.DataFrame:
        """Format and clean final annotation results.

        Args:
            df_list: List of annotation DataFrames.

        Returns:
            pd.DataFrame: Formatted merged results.
        """
        columns_to_drop = [
            "fragment_mz_delta",
            "fragment_fatty_acid_pos",
            "fragment_types"
        ]
        df_result_list = df_list.copy()
        df_result = pd.concat(
            df_result_list,
            axis=0,
            ignore_index=True
        )
        df_result = df_result.drop(columns=columns_to_drop)
        df_result = df_result.rename(columns={
                                    "precursor_mz": "mz_measured",
                                })
        df_result = self._compute_mz_error(data=df_result)
        df_result = self._move_after_columns(
                            df=df_result,
                            move_cols=[
                                "structure_annotated_name", 
                                "rt_score",
                                "mz_score",
                                "spec_score",
                                "total_score",
                                "fragments_ion_detection_level",
                                "fragments_spectrum_detection_level",
                            ],
                            target_col="primary_annotated_name"
                        )
        df_result = self._move_after_columns(
                            df=df_result,
                            move_cols=[
                                "smiles"
                            ],
                            target_col="mol_formula"
                        )
        
        return df_result
    
    def process_df_features_all(
        self,
        df: pd.DataFrame
    ) -> pd.DataFrame:
        """Format and clean all features (parent + fragments).

        Args:
            df: Merged features DataFrame.

        Returns:
            pd.DataFrame: Formatted full feature table.
        """
        df_features_all = df.copy()
        df_features_all = df_features_all.rename(columns={
                                    "precursor_mz": "mz_measured",
                                })
        df_features_all = self._fill_na_in_df_features_all(
                            df=df_features_all, 
                            na_col="mz_calculated"
                        )
        df_features_all["mz_calculated"] = df_features_all["mz_calculated"] + df_features_all["fragment_mz_delta"]
        df_features_all = self._fill_na_in_df_features_all(
                            df=df_features_all, 
                            na_col="adduct"
                        )
        df_features_all = self._fill_na_in_df_features_all(
                            df=df_features_all, 
                            na_col="mol_formula"
                        )
        # compute mz error (ppm)
        df_features_all = self._compute_mz_error(data=df_features_all) 
    
        df_features_all = self._move_after_columns(
                            df=df_features_all,
                            move_cols=[
                                "structure_annotated_name",
                                "rt_score",
                                "mz_score",
                                "spec_score",
                                "total_score",
                                "fragments_ion_detection_level",
                                "fragments_spectrum_detection_level",
                                "fragment_mz_delta",
                                "fragment_fatty_acid_pos",
                                "fragment_types",
                            ],
                            target_col="primary_annotated_name"
                        )
        df_features_all = self._move_after_columns(
                            df=df_features_all,
                            move_cols=[
                                "smiles",
                            ],
                            target_col="mol_formula"
                        )
        return df_features_all

    def _compute_mz_error(
        self, 
        data: pd.DataFrame
    ) -> pd.DataFrame:
        """Calculate mass error in ppm.

        Args:
            data: DataFrame with mz_measured and mz_calculated.

        Returns:
            pd.DataFrame: DataFrame with mz_error(ppm) column.
        """
        df = data.copy()
        df["mz_calculated"] = df["mz_calculated"].astype(float)
        ppm_series = abs(df["mz_measured"] - df["mz_calculated"]) / df["mz_calculated"] * 1e6 
        pos = df.columns.get_loc("mz_calculated") + 1
        df.insert(loc=pos, column="mz_error(ppm)", value=ppm_series)   
        
        return df
    
    def _move_after_columns(
        self,
        df: pd.DataFrame, 
        move_cols: list,
        target_col: str
    ) -> pd.DataFrame:
        """Reorder columns to place move_cols after target_col.

        Args:
            df: Input DataFrame.
            move_cols: Columns to move.
            target_col: Target position column.

        Returns:
            pd.DataFrame: Reordered DataFrame.
        """
        if isinstance(move_cols, str):
            move_cols = [move_cols]
        pos = df.columns.get_loc(target_col) + 1
        others = [c for c in df.columns if c not in move_cols]
        new_order = others[:pos] + move_cols + others[pos:]

        return df[new_order]
    
    def _fill_na_in_df_features_all(
        self,
        df: pd.DataFrame,
        na_col: str = "mz_calculated"
    ) -> pd.DataFrame:
        """Fill missing values using structure_annotated_name mapping.

        Args:
            df: Features DataFrame.
            na_col: Column to fill missing values.

        Returns:
            pd.DataFrame: DataFrame with filled values.
        """
        df_features_all = df.copy()
        fill_map = (
            df_features_all[df_features_all[na_col] != ""]
            .groupby("structure_annotated_name")[na_col]
            .first()
        )
        mask = df_features_all[na_col] == ""
        df_features_all.loc[mask, na_col] = (
            df_features_all.loc[mask, "structure_annotated_name"].map(fill_map)
        )
        return df_features_all
    
    def get_features_all(self) -> pd.DataFrame:
        """Return the full DataFrame of all annotated features (parent + fragments).

        Returns:
            pd.DataFrame: All annotated features.
        """
        return self.features_all
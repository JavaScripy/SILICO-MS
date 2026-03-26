
"""Core processing utilities for lipidomics LC-MS/MS feature analysis.

This module provides functions and classes for:
1. Splitting MS1 features into annotated candidates and unannotated auxiliaries
2. Mapping MS2 spectra to feature IDs for fast lookup
3. Checking lipid double bond status from annotation names
4. Spectrum-feature processing with RT/mz filtering and similarity scoring
5. Calculating multi-modal similarity (RT + m/z + MS2 spectrum)
6. Searching similar features using retention time and mass tolerance

Classes:
    SpectrumFeatureProcesser: Main processor for feature filtering and similarity calculation.

Functions:
    split_features: Separates candidate and auxiliary LC-MS features.
    ms2_spectra_to_dict: Creates a feature ID -> Spectrum lookup dictionary.
    get_spec_from_feature: Retrieves MS2 spectrum by feature ID.
    is_lipid_has_double_bond: Checks if lipid annotation indicates zero double bonds.
"""

from typing import List, Tuple, Dict
import re

import pandas as pd
from matchms.similarity import (
    BaseSimilarity, CosineGreedy,
    ModifiedCosine, NeutralLossesCosine,
    CosineHungarian
)
from matchms import Spectrum


def split_features(
    ms1_peak_table: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Splits MS1 peak table into candidate and auxiliary features.

    Candidate features have non-empty primary annotations.
    Auxiliary features have no primary annotations.

    Args:
        ms1_peak_table: DataFrame containing all LC-MS features.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]:
            DataFrames for candidate features and auxiliary features.
    """
    df_candidate_features = ms1_peak_table[ms1_peak_table["primary_annotated_name"] != ""]
    df_auxiliary_features = ms1_peak_table[ms1_peak_table["primary_annotated_name"] == ""]
    
    return df_candidate_features, df_auxiliary_features


def ms2_spectra_to_dict(
    ms2_spectra: List[Spectrum]
) -> Dict[str, Spectrum]:
    """Converts list of MS2 spectra into a feature ID lookup dictionary.

    Args:
        ms2_spectra: List of matchms Spectrum objects with feature_id metadata.

    Returns:
        Dict[str, Spectrum]: Dictionary mapping feature_id to Spectrum.
    """
    return {
        spec.get("feature_id"): spec
        for spec in ms2_spectra
    }


def get_spec_from_feature(
    feature: pd.Series,
    spectra_dict: dict
) -> Spectrum:
    """Retrieves the corresponding MS/MS spectrum for a given feature.

    Looks up the spectrum using the feature_id stored in the feature Series.

    Args:
        feature: Pandas Series containing a single LC-MS feature with feature_id.
        spectra_dict: Dictionary mapping feature_id strings to Spectrum objects.

    Returns:
        Spectrum: Matching MS/MS Spectrum if found; None otherwise.
    """
    feature_id = feature.get("feature_id")
    return spectra_dict.get(feature_id, None)


def is_lipid_has_double_bond(primary_annotated_name: str) -> bool:
    """Checks if lipid annotation indicates zero double bonds.

    Extracts numbers after colon (e.g., 18:0 → 0 double bonds)
    and returns True if all values are zero.

    Args:
        primary_annotated_name: Lipid annotation string.

    Returns:
        bool: True if lipid has no double bonds, False otherwise.
    """
    double_bond_list = re.findall(r':(.)', primary_annotated_name)
    double_bond_list = [
            int(double_bond) 
            for double_bond in double_bond_list
    ]
    return all(x == 0 for x in double_bond_list)


class SpectrumFeatureProcesser:
    """Handles LC-MS feature filtering, scoring, and similarity calculations.

    Provides configurable tolerance for retention time and m/z,
    supports multiple spectrum similarity algorithms,
    and computes combined similarity scores for lipid fragment matching.

    Attributes:
        rt_tol: Retention time tolerance threshold.
        rt_tol_mode: "absolute" or "relative" tolerance mode.
        mz_tol: Precursor m/z tolerance threshold.
        mz_tol_mode: "Da" or "ppm" tolerance unit.
        rt_weight: Weight for RT in total similarity score.
        precursor_mz_weight: Weight for m/z in total similarity score.
        spec_weight: Weight for MS2 spectrum in total similarity score.
        feature_sim_rt_weight: Weight for RT in fast feature search.
        feature_sim_mz_weight: Weight for m/z in fast feature search.
        spectrum_similarity: Configured matchms similarity metric.
    """
    def __init__(
        self,
        rt_tol: float = 0.5,
        rt_tol_mode: str = "absolute",
        mz_tol: float = 1.0,
        mz_tol_mode: str = "Da",
        rt_weight: float = 0.25,
        precursor_mz_weight: float = 0.25,
        spec_weight: float = 0.5,
        feature_sim_rt_weight: float = 1.0,
        feature_sim_mz_weight: float = 0.0,
        ms2_spectrum_similarity_type: str = "None"
    ):
        """Initializes feature processor with tolerance and weighting parameters.

        Args:
            rt_tol: Retention time tolerance.
            rt_tol_mode: "absolute" or "relative".
            mz_tol: Precursor m/z tolerance.
            mz_tol_mode: "Da" or "ppm".
            rt_weight: Weight of RT in total score.
            precursor_mz_weight: Weight of m/z in total score.
            spec_weight: Weight of MS2 spectrum in total score.
            feature_sim_rt_weight: RT weight for fast similar feature search.
            feature_sim_mz_weight: m/z weight for fast similar feature search.
            ms2_spectrum_similarity_type: Spectrum similarity method name.
        """
        self.rt_tol = rt_tol
        self.rt_tol_mode = rt_tol_mode
        self.mz_tol = mz_tol
        self.mz_tol_mode = mz_tol_mode
        self.rt_weight = rt_weight
        self.precursor_mz_weight = precursor_mz_weight
        self.spec_weight = spec_weight
        
        self.feature_sim_rt_weight = feature_sim_rt_weight
        self.feature_sim_mz_weight = feature_sim_mz_weight

        self.spectrum_similarity = self._similarity_type(
                                    similarity_type=ms2_spectrum_similarity_type,
                                    mz_tol=mz_tol
                                )

    def _similarity_type(
        self,
        similarity_type: str = "None",
        mz_tol: float = 0.1
    ) -> BaseSimilarity:
        """Initializes spectrum similarity metric from configuration string.

        Supports multiple cosine-based algorithms with preset parameters.

        Args:
            similarity_type: Name of similarity method.
            mz_tol: m/z tolerance for peak matching.

        Returns:
            BaseSimilarity: Initialized matchms similarity object.

        Raises:
            ValueError: If unknown similarity type is provided.
        """
        match similarity_type:
            case "NIST-LC":
                spectrum_similarity = CosineGreedy(
                                        tolerance=mz_tol,
                                        mz_power=1.3,
                                        intensity_power=0.53
                                    )
            case "NIST-GC":
                spectrum_similarity = CosineGreedy(
                                        tolerance=mz_tol,
                                        mz_power=3.0,
                                        intensity_power=0.6
                                    )
            case "SQRT":
                spectrum_similarity = CosineGreedy(
                                        tolerance=mz_tol,
                                        mz_power=0.0,
                                        intensity_power=0.5
                                    )
            case "MassBank":
                spectrum_similarity = CosineGreedy(
                                        tolerance=mz_tol,
                                        mz_power=2.0,
                                        intensity_power=0.5
                                    )
            case "None":
                spectrum_similarity = CosineGreedy(
                                        tolerance=mz_tol,
                                        mz_power=0.0,
                                        intensity_power=1.0
                                    )
            case "ModifiedCosine":
                spectrum_similarity = ModifiedCosine(
                                        tolerance=mz_tol,
                                        mz_power=1.3,
                                        intensity_power=0.53
                                    )
            case "NeutralLossesCosine":
                spectrum_similarity = NeutralLossesCosine(
                                        tolerance=mz_tol,
                                        mz_power=1.3,
                                        intensity_power=0.53
                                    )
            case "CosineHungarian":
                spectrum_similarity = CosineHungarian(
                                        tolerance=mz_tol,
                                        mz_power=1.3,
                                        intensity_power=0.53
                                    )
            case _:
                raise ValueError(f"Unmatched spectrum similarity method: {similarity_type}")

        return spectrum_similarity

    def filter_by_rt(
        self,
        df_features: pd.DataFrame,
        rt_theory: float
    ) -> pd.DataFrame:
        """Filters features by retention time tolerance.

        Args:
            df_features: DataFrame of LC-MS features.
            rt_theory: Target retention time to filter around.

        Returns:
            pd.DataFrame: Filtered features within RT tolerance.

        Raises:
            ValueError: If invalid rt_tol_mode is used.
        """
        rt1 = df_features.get("rt")
        rt2 = rt_theory

        match self.rt_tol_mode:
            case "relative":
                rt_diff = abs(rt1 - rt2) * 100 / rt2
            case "absolute":
                rt_diff = abs(rt1 - rt2)
            case _:
                raise ValueError(f"rt_tol_mode must be 'relative' or 'absolute', got {self.rt_tol_mode!r}")

        return df_features[rt_diff < self.rt_tol]

    def filter_by_precursor_mz(
        self,
        df_features: pd.DataFrame,
        precursor_mz_theory: float,
    ) -> pd.DataFrame:
        """Filters features by precursor m/z tolerance.

        Args:
            df_features: DataFrame of LC-MS features.
            precursor_mz_theory: Target precursor m/z to filter around.

        Returns:
            pd.DataFrame: Filtered features within m/z tolerance.

        Raises:
            ValueError: If invalid mz_tol_mode is used.
        """
        mz1 = df_features.get("precursor_mz")
        mz2 = precursor_mz_theory

        match self.mz_tol_mode:
            case "ppm":
                mz_diff = abs(mz2 - mz1) * 10^6 / mz2
            case "Da":
                mz_diff = abs(mz2 - mz1)
            case _:
                raise ValueError(
                    "Invalid `mz_tol_mode`. " \
                    "Only `ppm` and `Da` are valid for `mz_tol_mode`!"
                )

        return df_features[mz_diff < self.mz_tol]

    def feature_similarity(
        self,
        feature1: pd.Series,
        feature2: pd.Series,
        mz_delta: float,
        ms2_spectra_dict: dict
    ) -> pd.Series:
        """Calculates multi-modal similarity between two features.

        Computes RT, m/z, and MS2 spectrum similarities,
        then returns weighted total score.

        Args:
            feature1: First LC-MS feature.
            feature2: Second LC-MS feature.
            mz_delta: Mass delta for m/z comparison.
            ms2_spectra_dict: Spectrum lookup dictionary.

        Returns:
            pd.Series: Individual scores and total similarity score.
        """
        rt_sim = self._rt_similarity(
                    feature1=feature1,
                    feature2=feature2
                )
        mz_sim = self._precursor_mz_similarity(
                    feature1=feature1,
                    feature2=feature2,
                    mz_delta=mz_delta
                )
        spec_sim = self._spectra_similarity(
                    feature1=feature1,
                    feature2=feature2,
                    ms2_spectra_dict=ms2_spectra_dict,
                )
        feature_sim = rt_sim * self.rt_weight + \
                        mz_sim * self.precursor_mz_weight + \
                        spec_sim * self.spec_weight

        return pd.Series({
            "rt_score": rt_sim,
            "mz_score": mz_sim,
            "spec_score": spec_sim,
            "total_score": feature_sim,
        })
        
    def _rt_similarity(
        self,
        feature1: pd.Series,
        feature2: pd.Series
    ) -> float:
        """Calculates normalized retention time similarity.

        Returns 1.0 for perfect match, decreases linearly to 0 at tolerance limit.

        Args:
            feature1: First feature.
            feature2: Second feature.

        Returns:
            float: RT similarity score between 0 and 1.
        """
        rt1 = feature1.get("rt")
        rt2 = feature2.get("rt")

        match self.rt_tol_mode:
            case "relative":
                rt_diff = abs(rt1 - rt2) * 100 / rt2
            case "absolute":
                rt_diff = abs(rt1 - rt2)
            case _:
                raise ValueError(
                    f"rt_tol_mode must be 'relative' or 'absolute', got {self.rt_tol_mode}"
                )

        if rt_diff < self.rt_tol:
            rt_sim = 1 - abs(rt1 - rt2) / self.rt_tol
        else:
            rt_sim = 0.0
    
        return rt_sim

    def _precursor_mz_similarity(
        self,
        feature1: pd.Series,
        feature2: pd.Series,
        mz_delta: float,
    ) -> float:
        """Calculates normalized precursor m/z similarity with mass delta.

        Args:
            feature1: First feature.
            feature2: Second feature.
            mz_delta: Mass shift applied to feature2 m/z.

        Returns:
            float: m/z similarity score between 0 and 1.
        """
        mz1 = feature1.get("precursor_mz")
        mz2 = feature2.get("precursor_mz") + mz_delta

        match self.mz_tol_mode:
            case "ppm":
                mz_diff = abs(mz2 - mz1) * 10^6 / mz2
            case "Da":
                mz_diff = abs(mz2 - mz1)
            case _:
                raise ValueError(
                    "Invalid `mz_tol_mode`. " \
                    "Only `ppm` and `Da` are valid for `mz_tol_mode`!"
                )

        if mz_diff < self.mz_tol:
            mz_sim = 1 - mz_diff / self.mz_tol
        else:
            mz_sim = 0.0
    
        return mz_sim

    def _spectra_similarity(
        self,
        feature1: pd.Series,
        feature2: pd.Series,
        ms2_spectra_dict: dict
    ) -> float:
        """Calculates MS2 spectrum similarity between two features.

        Returns 0 if either spectrum is missing.

        Args:
            feature1: First feature.
            feature2: Second feature.
            ms2_spectra_dict: Lookup dict for MS2 spectra.

        Returns:
            float: Spectrum similarity score between 0 and 1.
        """
        spec1 = get_spec_from_feature(
                    feature=feature1, 
                    spectra_dict=ms2_spectra_dict
                )
        spec2 = get_spec_from_feature(
                    feature=feature2, 
                    spectra_dict=ms2_spectra_dict
                )

        if spec1 is None or spec2 is None:
            spec_sim = 0.0
        else:
            score = self.spectrum_similarity.pair(spec1, spec2)
            spec_sim =  float(score["score"])

        return spec_sim

    def search_similar_feature(
        self,
        feature: pd.Series,
        df_features_blank: pd.DataFrame,
    ) -> pd.DataFrame:
        """Searches for the most similar auxiliary feature.

        Filters by RT and m/z, computes similarity score,
        returns top-1 similar feature with annotation propagation.

        Args:
            feature: Query candidate feature.
            df_features_blank: Auxiliary (unannotated) features to search.

        Returns:
            pd.DataFrame: Top similar feature with propagated annotations.
        """
        rt, mz = feature.get("rt"), feature.get("precursor_mz")
        df = self.filter_by_rt(
                df_features=df_features_blank, 
                rt_theory=rt
            )
        if df.empty:
            return pd.DataFrame()
        
        df = self.filter_by_precursor_mz(
                df_features=df,
                precursor_mz_theory=mz
            )
        if df.empty:
            return pd.DataFrame()
        
        df["sim_score"] = self.feature_sim_rt_weight * (1 - abs(df["rt"] - rt) / rt)  + \
                          self.feature_sim_mz_weight * (1 - abs(df["precursor_mz"] - mz) / mz)
        df_similar_feature = df.nlargest(n=1, columns=["sim_score"])
        
        df_similar_feature = df_similar_feature.drop(columns=["sim_score"]) 
        df_similar_feature["fragment_mz_delta"] = feature.get("fragment_mz_delta", 0)
        df_similar_feature["primary_annotated_name"] = feature.get("primary_annotated_name", "")
        df_similar_feature["structure_annotated_name"] = feature.get("structure_annotated_name", "")
        df_similar_feature["fragment_fatty_acid_pos"] = feature.get("fragment_fatty_acid_pos", 0)
        df_similar_feature["fragment_types"] = feature.get("fragment_types", "")
        df_similar_feature["score"] = feature.get("score", 1.0)
        
        return df_similar_feature
    



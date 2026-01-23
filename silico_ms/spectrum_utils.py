
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
    df_candidate_features = ms1_peak_table[ms1_peak_table["primary_annotated_name"] != ""]
    df_auxiliary_features = ms1_peak_table[ms1_peak_table["primary_annotated_name"] == ""]
    
    return df_candidate_features, df_auxiliary_features

def ms2_spectra_to_dict(ms2_spectra: List[Spectrum]) -> Dict[str, Spectrum]:
    return {
            spec.get("feature_id"): spec
            for spec in ms2_spectra
    }

def get_spec_from_feature(
        feature: pd.Series,
        spectra_dict: dict
    ) -> Spectrum:
    """Get spectrum by feature ID.

    Parameters:
        feature:
            LC-MS feature.
        spectra:
            All of MS/MS spectra.

    Returns:
        MS/MS Spectrum of the feature.
    """
    feature_id = feature.get("feature_id")
    
    return spectra_dict.get(feature_id, None)

def is_lipid_has_double_bond(primary_annotated_name: str) -> bool:
    """"""
    double_bond_list = re.findall(r':(.)', primary_annotated_name)
    double_bond_list = [
            int(double_bond) 
            for double_bond in double_bond_list
    ]
    return all(x == 0 for x in double_bond_list)


class SpectrumFeatureProcesser:
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
        """"""

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
        """Calculate similarity between two spectrum feature."""
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

        #return feature_sim
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
        """Calculate similarity of retention time (RT)."""
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
        """Calculate similarity of precursor m/z."""
        
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
        """Calculate similarity of MS/MS spctrum."""

        spec1 = get_spec_from_feature(
                    feature=feature1, 
                    spectra_dict=ms2_spectra_dict
                )
        spec2 = get_spec_from_feature(
                    feature=feature2, 
                    spectra_dict=ms2_spectra_dict
                )
        
        #print("spec1: ", spec1)
        #print("spec2: ", spec2)

        if spec1 is None or spec2 is None:
            spec_sim = 0.0
        else:
            score = self.spectrum_similarity.pair(spec1, spec2)
            spec_sim =  float(score["score"])
        
        #if spec_sim > 0.4:
        #    print("spec_sim:", spec_sim)

        return spec_sim

    def search_similar_feature(
        self,
        feature: pd.Series,
        df_features_blank: pd.DataFrame,
    ) -> pd.DataFrame:
        """
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
    



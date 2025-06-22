from typing import List, Union


import numpy as np 
import matchms
from matchms import calculate_scores, Spectrum, Scores
from matchms.similarity import CosineGreedy, ModifiedCosine, NeutralLossesCosine, CosineHungarian
from matchms.similarity import PrecursorMzMatch, MetadataMatch






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

# need edit
def calculate_similarity_scores(queries: Union[List[Spectrum], Spectrum],
                                references: List[Spectrum],
                                score_type: str = "None",
                                mz_tol: float = 0.1,
                                *args,
                                **kwargs) -> Scores:
    similarity_score = get_similarity_type(score_type, mz_tol)
    if not isinstance(queries, list):
        queries = [queries]
    try:
        scores = calculate_scores(references, queries, similarity_score)
    except:
        scores = None
    return scores


def get_similarity_vector(scores: Scores) -> np.ndarray:
    scores_vector = scores.to_array().squeeze()
    similarity = np.array([t[0] for t in scores_vector])
    return similarity


def get_matched_peaks_vector(scores: Scores) -> np.ndarray:
    scores_vector = scores.to_array().squeeze()
    matched_peaks = np.array([t[1] for t in scores_vector])
    return matched_peaks


def get_idx_filter_by_precursor_mz(queries: Union[List[Spectrum], Spectrum],
                                   references: List[Spectrum],
                                   mz_tol: float = 0.1,
                                   tolerance_mode = "Dalton") -> np.ndarray:
    """
    match tolerance_mode:
        case "Dalton":
            print("Unit of mass tolerance is Dalton.")
        case "ppm":
            print("Unit of mass tolerance is ppm.")        
        case _:
            raise ValueError("Invalid mode, mode of tolerance must be 'Dalton' or 'ppm'")
    """
    
    if not isinstance(queries, list):
        queries = [queries]
    similarity_score = PrecursorMzMatch(mz_tol, tolerance_mode)
    scores = calculate_scores(references, queries, similarity_score)
    scores_vector = scores.scores.to_array().squeeze()
    return np.argwhere(scores_vector).squeeze()


def get_idx_filter_by_metadata(queries: Union[List[Spectrum], Spectrum],
                               references: List[Spectrum],
                               field: str = "ionmode",
                               matching_type: str = "euqal_match",
                               tolerance: float = 0.1) -> np.ndarray:
    if not isinstance(queries, list):
        queries = [queries]
    similarity_score = MetadataMatch(field, matching_type, tolerance)
    scores = calculate_scores(references, queries, similarity_score)
    scores_vector = scores.scores.to_array().squeeze()
    return np.argwhere(scores_vector).squeeze()
    

def get_idx_filter_by_similarity_score(socres: Scores,
                                       score_thred: float = 0.6):
    pass

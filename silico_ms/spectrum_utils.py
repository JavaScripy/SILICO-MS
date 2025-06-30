
from typing import Tuple, Optional, Dict
from matchms import Spectrum

from matchms.similarity import CosineGreedy, ModifiedCosine, NeutralLossesCosine, CosineHungarian


def get_rt_similarity(
        feature1: dict,
        feature2: dict,
        tolerance: float = 0.5,
        mode: str = "absolute"
    ):
    """Get similarity of retention time (RT).

    Parameters:
        feature1:
            Feature of query.
        feature2:
            Feature of reference.
        tolerance:
            Tolerance of RT.
        mode:
            Type of tolerance.

    Returns:
        rt_score:
            Score of RT.
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


def get_precursor_mz_similarity(
        feature1: dict,
        feature2: dict,
        tolerance: float = 5.0,
        mode: str = "ppm"
    ) -> float:
    """Get similarity of m/z.

    Parameters:
        feature1:
            Feature of query.
        feature2:
            Feature of reference.
        tolerance:
            Tolerance of m/z.
        mode:
            Type of tolerance.

    Returns:
        mz_score:
            Score of m/z.
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


def get_ms2_spec_similarity(
        feature1: dict,
        feature2: dict,
        spectra: Dict[str, Spectrum],
        mz_tol: float = 0.5,
        score_type: str = "None"
    ) -> Tuple[float, int]:
    """Get similarity of MS/MS spectrum.

    Parameters:
        feature1:
            Feature of query.
        feature2:
            Feature of reference.
        spectra:
            All of MS/MS spectra.
        tolerance:
            Tolerance of m/z.
        mode:
            Type of tolerance.

    Returns:
        spec_score:
            Score of spectrum similarity score.
        spec_match:
            Number of matched peaks.
    """
    similarity_score = get_similarity_type(score_type, mz_tol)
    spec1 = get_spec_from_feature(feature1, spectra)
    spec2 = get_spec_from_feature(feature2, spectra)

    if spec1 is None or spec2 is None:
        return 0.0
    else:
        score = similarity_score.pair(spec1, spec2)
        spec_score =  float(score["score"])
        spec_match = int(score["matches"])
        return spec_score, spec_match


def get_similarity_type(
        score_type: str = "None",
        mz_tol: float = 0.1
    ):
    """Get the type of MS/MS similarity type.

    Parameters:
        score_type:
            Type of similarity score.
        mz_tol:
            Tolerance of m/z.

    Returns:
        similarity_score:
            Type of similarity score.
    """
    match score_type:
        case "NIST-LC":
            similarity_score = CosineGreedy(tolerance=mz_tol, mz_power=1.3, intensity_power=0.53)
        case "NIST-GC":
            similarity_score = CosineGreedy(tolerance=mz_tol, mz_power=3.0, intensity_power=0.6)
        case "SQRT":
            similarity_score = CosineGreedy(tolerance=mz_tol, mz_power=0.0, intensity_power=0.5)
        case "MassBank":
            similarity_score = CosineGreedy(tolerance=mz_tol, mz_power=2.0, intensity_power=0.5)
        case "None":
            similarity_score = CosineGreedy(tolerance=mz_tol, mz_power=0.0, intensity_power=1.0)
        case "ModifiedCosine":
            similarity_score = ModifiedCosine(tolerance=mz_tol, mz_power=1.3, intensity_power=0.53)
        case "NeutralLossesCosine":
            similarity_score = NeutralLossesCosine(tolerance=mz_tol, mz_power=1.3, intensity_power=0.53)
        case "CosineHungarian":
            similarity_score = CosineHungarian(tolerance=mz_tol, mz_power=1.3, intensity_power=0.53)
        case _:
            print("Not matched method!")

    return similarity_score


def get_spec_from_feature(
        feature: dict,
        spectra: dict
    ) -> Optional[Spectrum]:
    """Get spectrum by feature ID.

    Parameters:
        feature:
            LC-MS feature.
        spectra:
            All of MS/MS spectra.

    Returns:
        MS/MS Spectrum of the feature.
    """
    feature_id = feature["id"]
    return spectra.get(feature_id, None)


def filter_by_rt(
        feature: dict,
        rt_theory: float,
        rt_tol: float = None,
        rt_mode: str = "absolute"
    ) -> Optional[dict]:
    """Filter feature by retention time (RT).

    Parameters:
        feature:
            Feature of query.
        rt_theory:
            Theoretical RT.
        rt_tol:
            Tolerance of RT.
        rt_mode:
            Mode of tolerance of RT.

    Returns:
        feature:
            Feature of output.
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


def filter_by_mz(
        feature: dict,
        mz_theory: float,
        mz_tol: float = None,
        mz_mode: str = "Da"
    ) -> dict:
    """Filter feature by m/z.

    Parameters:
        feature:
            Feature of query.
        mz_theory:
            Theoretical m/z.
        mz_tol:
            Tolerance of mz.
        mz_mode:
            Mode of tolerance of m/z.

    Returns:
        feature:
            Feature of output.
    """
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

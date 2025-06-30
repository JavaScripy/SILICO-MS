
from matchms.similarity import CosineGreedy, ModifiedCosine, NeutralLossesCosine, CosineHungarian






def get_rt_similarity(feature1: dict, 
                      feature2: dict, 
                      tolerance: float = 0.5, 
                      mode: str = "absolute"):
    """edited from MetDNA.
    Args:
        feature1: query
        feature2: reference
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


def get_precursor_mz_similarity(feature1: dict, 
                                feature2: dict, 
                                tolerance = 5, 
                                mode = "ppm"):
    """edited from MetDNA.
    Args:
        feature1: query
        feature2: reference
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


def get_ms2_spec_similarity(feature1: dict, 
                            feature2: dict, 
                            spectra: dict, 
                            mz_tol: float = 0.5, 
                            score_type: str = "None"):
    similarity_score = get_similarity_type(score_type, mz_tol)
    spec1 = get_spec_from_feature(feature1, spectra)
    spec2 = get_spec_from_feature(feature2, spectra)
    if spec1 is None or spec2 is None:
        return 0.0
    else:
        score = similarity_score.pair(spec1, spec2)
        return score['score'], score['matches']


def get_similarity_type(score_type: str = "None",
                        mz_tol: float = 0.1):
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


def get_spec_from_feature(feature: dict, 
                          spectra: dict):
    feature_id = feature["id"]
    return spectra.get(feature_id, None)


def filter_by_rt(feature: dict, 
                 rt_theory: float, 
                 rt_tol: float = None, 
                 rt_mode: str = "absolute"):
    """filter by Retention Time
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
        

def filter_by_mz(feature: dict, 
                 mz_theory: float, 
                 mz_tol: float = None, 
                 mz_mode: str = "Da"):
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
    
from typing import List
from itertools import combinations_with_replacement, product, chain
import re


import pandas as pd




def generate_lipids_batch(
    fatty_acids: List[str],
    lipid_classes: List[dict],
    ceramide_bases: List[str] = None,
    sterol_bases: List[str] = None
) -> List[str]:
    """"""
    lipids = [
        generate_lipid_single(
            fatty_acids=fatty_acids,
            lipid_class=lipid_class,
            ceramide_bases=ceramide_bases,
            sterol_bases=sterol_bases
        )
        for lipid_class in lipid_classes
    ]
    lipids = list(chain.from_iterable(lipids))
    lipids = [
        lipid 
        for lipid in lipids 
        if lipid is not None
    ]
    return lipids
    


def generate_lipid_single(
    fatty_acids: List[str],
    lipid_class: dict,
    ceramide_bases: List[str] = None,
    sterol_bases: List[str] = None,
    #num_fatty_acid: int
):
    lipid_class_name = lipid_class.get("lipid_class")
    num_acyl_chain = lipid_class.get("num_acyl_chain")
    num_acyl_chain = int(num_acyl_chain)
    is_ceramide_based = lipid_class.get("is_ceramide_based")
    is_sterol_based = lipid_class.get("is_sterol_based")
    
    if is_sterol_based:
        lipids = _generate_sterol_base_liipids(
            fatty_acids=fatty_acids,
            lipid_class=lipid_class,
            sterol_bases=sterol_bases
        )
    elif is_ceramide_based:
        lipids = _generate_ceramide_base_lipids(
                    ceramide_bases=ceramide_bases,
                    fatty_acids=fatty_acids,
                    lipid_class_name=lipid_class_name
                )    
    else:
        not_ether_lipids = _generate_normal_lipids(
                    lipid_class_name=lipid_class_name,
                    num_acyl_chain=num_acyl_chain,
                    fatty_acids=fatty_acids,
                )
        ether_lipids = _generate_ether_lipids(
                    lipid_class_name=lipid_class_name,
                    num_acyl_chain=num_acyl_chain,
                    fatty_acids=fatty_acids
                )
        plasmanyl_lipids = _generate_plasmanyl_ether_lipids(
                    lipid_class_name=lipid_class_name,
                    num_acyl_chain=num_acyl_chain,
                    fatty_acids=fatty_acids
                )
        lipids = not_ether_lipids + ether_lipids + plasmanyl_lipids
        
    return lipids
    

def _generate_normal_lipids(
    lipid_class_name: str,
    num_acyl_chain: int,
    fatty_acids: List[str]
) -> List[str]:
    """
    """
    total_fatty_acids = list(combinations_with_replacement(fatty_acids, num_acyl_chain))
    lipids = [
        f"{lipid_class_name} {'_'.join(total_fatty_acid)}"
        for total_fatty_acid in total_fatty_acids
    ]
    return lipids

def _generate_ether_lipids(
    lipid_class_name: str,
    num_acyl_chain: int,
    fatty_acids: List[str]
) -> List[str]:
    not_ether_fatty_acids = list(combinations_with_replacement(fatty_acids, num_acyl_chain-1))
    ether_fatty_acids = [
        f"O-{fatty_acid}" 
        for fatty_acid in fatty_acids
    ]
    total_fatty_acids =  [
        (ether_fa, ) + not_ether_fa
        for ether_fa, not_ether_fa
        in product(ether_fatty_acids, not_ether_fatty_acids)
    ]
    lipids = [
        f"{lipid_class_name} {'_'.join(total_fatty_acid)}"
        for total_fatty_acid in total_fatty_acids
    ]
    return lipids


def _generate_plasmanyl_ether_lipids(
    lipid_class_name: str,
    num_acyl_chain: int,
    fatty_acids: List[str]
) -> List[str]:
    if lipid_class_name not in ["PC", "PE"]:
        return []

    not_ether_fatty_acids = list(combinations_with_replacement(fatty_acids, num_acyl_chain-1))
    ether_fatty_acids = [
        f"P-{fatty_acid}" 
        for fatty_acid in fatty_acids
    ]
    total_fatty_acids =  [
        (ether_fa, ) + not_ether_fa
        for ether_fa, not_ether_fa
        in product(ether_fatty_acids, not_ether_fatty_acids)
    ]
    lipids = [
        f"{lipid_class_name} {'_'.join(total_fatty_acid)}"
        for total_fatty_acid in total_fatty_acids
    ]
    return lipids


def _generate_ceramide_base_lipids(
    lipid_class_name: str,
    fatty_acids: List[str],
    ceramide_bases: List[str],
):
    """
    """
    total_fatty_acids = product(ceramide_bases, fatty_acids)
    lipids = [
        f"{lipid_class_name} {'_'.join(total_fatty_acid)}"
        for total_fatty_acid in total_fatty_acids
    ]
    return lipids

def _generate_sterol_base_liipids(
    fatty_acids: List[str],
    lipid_class_name: str,
    sterol_bases: List[str] = None,   
):
    if sterol_bases is None:
        return None
    
    total_fatty_acids = product(sterol_bases, fatty_acids)
    lipids = [
        f"{lipid_class_name} {'_'.join(total_fatty_acid)}"
        for total_fatty_acid in total_fatty_acids
    ]
    return lipids


def db_shorthand_to_chain_shorthand(
    db_shorthands: List[str] 
) -> List[str]:
    """
    """
    chain_shorthands = [
        _db_shorthand_to_chain_shorthand_base(db_shorthand=shorthand)
        for shorthand in db_shorthands
    ]
    return chain_shorthands


def _db_shorthand_to_chain_shorthand_base(
    db_shorthand: str
):
    DB_PATTERN = r"\([\d,ZE]+\)"
    chain_shorthand = re.sub(
        pattern=DB_PATTERN,
        repl="",
        string=db_shorthand
    )
    return chain_shorthand

    
if __name__ == "__main__":
    
    fatty_acid_file = "raw_data/all_cuatred_acyl_fatty_acid.csv"
    cer_base_file = "raw_data/all_cuatred_ceramide_base_fatty_acid.csv"
    lipid_class_file = "raw_data/all_cuatred_lipid_class.csv"
    out_file = "clean_data/lipidmaps_chain_shorthand.csv"
    
    df_fattty_acid = pd.read_csv(fatty_acid_file)
    df_lipid_class = pd.read_csv(lipid_class_file)
    df_ceramide_base = pd.read_csv(cer_base_file)
    fatty_acid_list = df_fattty_acid["fatty_acid"].tolist()
    fatty_acid_list = [
        fatty_acid.split(" ", 1)[-1]
        for fatty_acid in fatty_acid_list
    ]
    lipid_class_list = df_lipid_class.to_dict(orient="records")
    ceramide_base_list = df_ceramide_base["fatty_acid"].tolist()
    ceramide_base_list = [
        ceramide_base.split(" ", 1)[-1]
        for ceramide_base in ceramide_base_list
    ]
    
    
    db_shorthands = generate_lipids_batch(
                        fatty_acids=fatty_acid_list,
                        ceramide_bases=ceramide_base_list,
                        lipid_classes=lipid_class_list,
                    )
    chain_shorthands = db_shorthand_to_chain_shorthand(db_shorthands=db_shorthands)
    
    df = pd.DataFrame(
        data={
            "chain_shorthand": chain_shorthands,
            "db_shorthand": db_shorthands  
        })
    df.to_csv(out_file, index=None)
    
    
    
    
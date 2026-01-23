from typing import List
from itertools import combinations_with_replacement, product, chain
import re


from tqdm import tqdm
import pandas as pd


class LipidGenerator:    

    def __init__(
        self,
        df_fatty_acids: pd.DataFrame,
        df_lipid_classes: pd.DataFrame,
        df_ceramide_bases: pd.DataFrame
    ):
        self.fatty_acids = self._init_fatty_acids(df_fatty_acids)
        self.lipid_classes = self._init_lipid_class(df_lipid_classes)
        self.ceramide_bases = self._init_ceramide_bases(df_ceramide_bases)

    def _init_fatty_acids(
        self, 
        df_fatty_acids: pd.DataFrame
    ) -> List[str]:
        """"""
        fatty_acid_list = df_fatty_acids["fatty_acid"].astype(str).tolist()
        fatty_acid_list = [
            fatty_acid.split(" ", 1)[-1]
            for fatty_acid in fatty_acid_list
        ]
        fatty_acid_list = [fa.strip() for fa in fatty_acid_list]
        return fatty_acid_list
    
    def _init_lipid_class(
        self,
        df_lipid_classes: pd.DataFrame,
    ) -> List[dict]:
        """"""
        return df_lipid_classes.to_dict(orient="records")
    
    def _init_ceramide_bases(
        self,
        df_ceramide_bases: pd.DataFrame
    ) -> List[dict]:
        """"""
        ceramide_fa_list = df_ceramide_bases["fatty_acid"].astype(str).tolist()
        ceramide_fa_list = [
            ceramide_fa.split(" ", 1)[-1]
            for ceramide_fa in ceramide_fa_list
        ]
        ceramide_fa_list = [fa.strip() for fa in ceramide_fa_list]
        ceramide_only_list = df_ceramide_bases["ceramide_only"].tolist()
        ceramide_bases = [
            {"fatty_acid": fa, "ceramide_only": only} 
            for fa, only in zip(ceramide_fa_list, ceramide_only_list)
        ]
        return ceramide_bases
    
    def generate_df_lipids(self) -> pd.DataFrame:
        """Generate a DataFrame of lipids."""
        db_shorthands = self.generate_lipids_batch()
        chain_shorthands = self.db_shorthand_to_chain_shorthand(db_shorthands=db_shorthands)

        df_lipids = pd.DataFrame(
            data={
            "chain_shorthand": chain_shorthands,
            "db_shorthand": db_shorthands  
        })
        
        return df_lipids
        
    def generate_lipids_batch(self) -> List[str]:
        """Generate a list of lipids."""
        lipids = [
            self._generate_lipid_single(
                fatty_acids=self.fatty_acids,
                lipid_class=lipid_class,
                ceramide_bases=self.ceramide_bases
            )
            for lipid_class in tqdm(self.lipid_classes, desc="Generating lipids")
        ]
        lipids = list(chain.from_iterable(lipids))
        lipids = [
            lipid 
            for lipid in lipids 
            if lipid is not None
        ]
        return lipids
    
    def _generate_lipid_single(
        self,
        fatty_acids: List[str],
        lipid_class: dict,
        ceramide_bases: List[dict] = None
    ) -> List[str]:
        """Generate a list of lipids."""
        lipid_class_name = lipid_class.get("lipid_class")
        num_acyl_chain = lipid_class.get("num_acyl_chain")
        num_acyl_chain = int(num_acyl_chain)
        is_ceramide_based = lipid_class.get("is_ceramide_based")
    
        if is_ceramide_based:
            lipids = self._generate_ceramide_base_lipids(
                    ceramide_bases=ceramide_bases,
                    fatty_acids=fatty_acids,
                    lipid_class_name=lipid_class_name
                )    
        else:
            not_ether_lipids = self._generate_normal_lipids(
                    lipid_class_name=lipid_class_name,
                    num_acyl_chain=num_acyl_chain,
                    fatty_acids=fatty_acids,
                )
            ether_lipids = self._generate_alkyl_ether_lipids(
                    lipid_class_name=lipid_class_name,
                    num_acyl_chain=num_acyl_chain,
                    fatty_acids=fatty_acids
                )
            plasmanyl_lipids = self._generate_plasmanyl_ether_lipids(
                    lipid_class_name=lipid_class_name,
                    num_acyl_chain=num_acyl_chain,
                    fatty_acids=fatty_acids
                )
            lipids = not_ether_lipids + ether_lipids + plasmanyl_lipids
        
        return lipids
    
    def _generate_normal_lipids(
        self,
        lipid_class_name: str,
        num_acyl_chain: int,
        fatty_acids: List[str]
    ) -> List[str]:
        """Generate a list of normal lipids."""
        total_fatty_acids = list(combinations_with_replacement(fatty_acids, num_acyl_chain))
        lipids = [
            f"{lipid_class_name} {'_'.join(total_fatty_acid)}"
            for total_fatty_acid in total_fatty_acids
        ]
        return lipids
    
    def _generate_alkyl_ether_lipids(
        self,
        lipid_class_name: str,
        num_acyl_chain: int,
        fatty_acids: List[str]
    ) -> List[str]:
        """Generate a list of ether lipids."""
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
        self,
        lipid_class_name: str,
        num_acyl_chain: int,
        fatty_acids: List[str]
    ) -> List[str]:
        """Generate a list of plasmanyl ether lipids."""
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
        self,
        lipid_class_name: str,
        fatty_acids: List[str],
        ceramide_bases: List[dict]
    ) -> List[str]:
        """Generate a list of ceramide-based lipids."""
        if lipid_class_name in ["Cer"]:
            base_fatty_acids = [
                base.get("fatty_acid") 
                for base in ceramide_bases 
            ]
        else:
            base_fatty_acids = [
                base.get("fatty_acid") 
                for base in ceramide_bases 
                if not base.get("ceramide_only")
            ]
        
        total_fatty_acids = product(base_fatty_acids, fatty_acids)
        lipids = [
            f"{lipid_class_name} {'_'.join(total_fatty_acid)}"
            for total_fatty_acid in total_fatty_acids
        ]
        return lipids
    
    def db_shorthand_to_chain_shorthand(
        self,
        db_shorthands: List[str] 
    ) -> List[str]:
        """Convert database shorthands to chain shorthands."""
        chain_shorthands = [
            self._db_shorthand_to_chain_shorthand_base(db_shorthand=shorthand)
            for shorthand in db_shorthands
        ]
        return chain_shorthands

    def _db_shorthand_to_chain_shorthand_base(
        self,
        db_shorthand: str
    ):
        """Convert a database shorthand to a chain shorthand."""
        DB_PATTERN = r"\([\d,ZE]+\)"
        chain_shorthand = re.sub(
            pattern=DB_PATTERN,
            repl="",
            string=db_shorthand
        )
        return chain_shorthand


def extra_CL_from_lipidmaps(data: pd.DataFrame) -> pd.DataFrame:
    """Process LipidMaps data and generate lipid list."""
    DB_PATTERN = r"\d+:\d+(\(\d+[ZE](,\d+[ZE])*\))*"
    CL_PATTERN = fr"^CL\(1\'-\[{DB_PATTERN}/{DB_PATTERN}\],3\'-\[{DB_PATTERN}/{DB_PATTERN}\]\)"
    
    df = data.copy()
    mask = df["COMMON_NAME"].str.match(CL_PATTERN, na=False)
    
    df_CL = df.loc[mask].copy()
    df_CL["db_shorthand"] = df_CL["COMMON_NAME"]
    df_CL["db_shorthand"] = df_CL["db_shorthand"].str.replace(pat=r"\[rac\]", repl="", regex=True)
    df_CL["db_shorthand"] = df_CL["db_shorthand"].str.replace(pat=r"\(1\'-\[", repl=" ", regex=True)
    df_CL["db_shorthand"] = df_CL["db_shorthand"].str.replace(pat=r"\]\)", repl="", regex=True)
    df_CL["db_shorthand"] = df_CL["db_shorthand"].str.replace(pat=r"\]\,3\'-\[", repl="_", regex=True)
    df_CL["db_shorthand"] = df_CL["db_shorthand"].str.replace(pat=r"\/", repl="_", regex=True)
    df_CL["db_shorthand"] = df_CL["db_shorthand"].str.replace(pat=r"[ZE]", repl="", regex=True)
    df_CL["chain_shorthand"] = df_CL["db_shorthand"].str.replace(pat=r"\(\d+(,\d+)*\)", repl="", regex=True)
    
    df_lipidmaps_CL = df_CL[["chain_shorthand", "db_shorthand"]]
    
    return df_lipidmaps_CL


if __name__ == "__main__":
    """"""
    fatty_acid_file = "raw_data/all_cuatred_acyl_fatty_acid.csv"
    cer_base_file = "raw_data/all_cuatred_ceramide_base_fatty_acid.csv"
    lipid_class_file = "raw_data/all_cuatred_lipid_class.csv"
    lipidmaps_file = "raw_data/lipidmaps.tsv"
    out_file = "clean_data/lipidmaps_chain_shorthand.csv"
    
    df_fattty_acid = pd.read_csv(fatty_acid_file)
    df_lipid_class = pd.read_csv(lipid_class_file)
    df_ceramide_base = pd.read_csv(cer_base_file)
    lipid_generator = LipidGenerator(
                        df_fatty_acids=df_fattty_acid,
                        df_lipid_classes=df_lipid_class,
                        df_ceramide_bases=df_ceramide_base
                    )
    df_lipids = lipid_generator.generate_df_lipids()
    
    df_lipidmaps = pd.read_csv(lipidmaps_file, sep="\t")
    df_lipidmaps_CL = extra_CL_from_lipidmaps(data=df_lipidmaps)
    
    df_lipids_total = pd.concat([df_lipids, df_lipidmaps_CL], axis=0)
    df_lipids_total.to_csv(out_file, index=None)


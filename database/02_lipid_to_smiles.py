from typing import List, Dict, Tuple
import re


from tqdm import tqdm
import pandas as pd

tqdm.pandas(desc="Lipid to SMILES")

"""
SMILES_TEMPLATES = {
    "PA": "{sn2}OC(CO{sn1})COP(=O)(O)O",  # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)O"
    "PC": "{sn2}OC(CO{sn1})COP(=O)([O-])OCC[N+](C)(C)C", # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)([O-])OCC[N+](C)(C)C"
    "PE": "{sn2}OC(CO{sn1})COP(=O)(O)OCCN", # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCCN"
    "PI": "{sn2}OC(CO{sn1})COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O"
    "PS": "{sn2}OC(CO{sn1})COP(=O)(O)OCC(N)C(=O)O", # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCC(N)C(=O)O"
    "PG": "{sn2}OC(CO{sn1})COP(=O)(O)OCC(O)CO", # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCC(O)CO"
    "MG": "{sn1}OCC(O)CO", # f"{sn1_forward}OCC(O)CO"
    "DG": "{sn2}OC(CO)CO{sn1}", # f"{sn2_forward}OC(CO)CO{sn1_reverse}"
    "TG": "{sn2}OC(CO{sn3})CO{sn1}", # f"{sn2_forward}OC(CO{sn3_reverse})CO{sn1_reverse}"
    "Cer": "{sn1}{sn2}", # f"{sn1_forward}{sn2_reverse}" # fa chain specific parse
    "SM": "{sn1}C(O)C(COP(=O)([O-])OCC[N+](C)(C)C)N{sn2}", # f"{sn1_forward}C(O)C(COP(=O)([O-])OCC[N+](C)(C)C)N{sn2_reverse}"
    "BMP": "{sn2}OC(CO)COP(=O)(O)OCC(CO)O{sn1}", # f"{sn2_forward}OC(CO)COP(=O)(O)OCC(CO)O{sn1_reverse}"
    "CAR": "{sn1}OC(CC(=O)[O-])C[N+](C)(C)C", # f"{sn1_forward}OC(CC(=O)[O-])C[N+](C)(C)C"
    "LPA": "{sn1}OCC(O)COP(=O)(O)O", # f"{sn1_forward}OCC(O)COP(=O)(O)O"
    "LPC": "{sn1}OCC(O)COP(=O)([O-])OCC[N+](C)(C)C", # f"{sn1_forward}OCC(O)COP(=O)([O-])OCC[N+](C)(C)C"
    "LPE": "{sn1}OCC(O)COP(=O)(O)OCCN", # f"{sn1_forward}OCC(O)COP(=O)(O)OCCN"
    "LPI": "{sn1}OCC(O)COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", # f"{sn1_forward}OCC(O)COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O"
    "LPS": "{sn1}OCC(O)COP(=O)(O)OCC(N)C(=O)O", # f"{sn1_forward}OCC(O)COP(=O)(O)OCC(N)C(=O)O"
    "LPG": "{sn1}OCC(O)COP(=O)(O)OCC(O)CO", # f"{sn1_forward}OCC(O)COP(=O)(O)OCC(O)CO"
    "CerP": "{sn1}C(O)C(COP(=O)(O)O)N{sn2}", # f"{sn1_forward}C(O)C(COP(=O)(O)O)N{sn2_reverse}"
    "CE": "{sn1}OC1CCC2(C)C(=CCC3C2CCC2(C)C(C(C)CCCC(C)C)CCC32)C1", # f"{sn1_forward}OC1CCC2(C)C(=CCC3C2CCC2(C)C(C(C)CCCC(C)C)CCC32)C1"
    "CL": "{sn4}OCC(COP(=O)(O)OCC(O)COP(=O)(O)OCC(CO{sn1})O{sn2})O{sn3}", # f"{sn4_forward}OCC(COP(=O)(O)OCC(O)COP(=O)(O)OCC(CO{sn1_reverse})O{sn2_reverse})O{sn3_reverse}"
    # ...
}

# reverse: omega-position
# forward: carboxyl-position
FA_CHAIN_REVERSE = {
    "PA": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)O"
    "PC": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)([O-])OCC[N+](C)(C)C"
    "PE": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCCN"
    "PI": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O"
    "PS": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCC(N)C(=O)O"
    "PG": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCC(O)CO"
    "MG": [False],          # f"{sn1_forward}OCC(O)CO"
    "DG": [True, False],    # f"{sn2_forward}OC(CO)CO{sn1_reverse}"
    "TG": [True, False, True],  # f"{sn2_forward}OC(CO{sn3_reverse})CO{sn1_reverse}"
    "Cer": [False, True],   # f"{sn1_forward}{sn2_reverse}" # fa chain specific parse
    "SM": [False, True],    # f"{sn1_forward}C(O)C(COP(=O)([O-])OCC[N+](C)(C)C)N{sn2_reverse}"
    "BMP": [True, False],   # f"{sn2_forward}OC(CO)COP(=O)(O)OCC(CO)O{sn1_reverse}"
    "CAR": [False],         # f"{sn1_forward}OC(CC(=O)[O-])C[N+](C)(C)C"
    "LPA": [False],         # f"{sn1_forward}OCC(O)COP(=O)(O)O"
    "LPC": [False],         # f"{sn1_forward}OCC(O)COP(=O)([O-])OCC[N+](C)(C)C"
    "LPE": [False],         # f"{sn1_forward}OCC(O)COP(=O)(O)OCCN"
    "LPI": [False],         # f"{sn1_forward}OCC(O)COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O"
    "LPS": [False],         # f"{sn1_forward}OCC(O)COP(=O)(O)OCC(N)C(=O)O"
    "LPG": [False],         # f"{sn1_forward}OCC(O)COP(=O)(O)OCC(O)CO"
    "CerP": [False, True],  # f"{sn1_forward}C(O)C(COP(=O)(O)O)N{sn2_reverse}" 
    "CE": [False],          # f"{sn1_forward}OC1CCC2(C)C(=CCC3C2CCC2(C)C(C(C)CCCC(C)C)CCC32)C1"
    "CL": [True, True, True, False],    # f"{sn4_forward}OCC(COP(=O)(O)OCC(O)COP(=O)(O)OCC(CO{sn1_reverse})O{sn2_reverse})O{sn3_reverse}"
    # ...
} 
"""


class LipidSmilesGenerator:
    SMILES_TEMPLATES = {
        "PA": "{sn2}OC(CO{sn1})COP(=O)(O)O",  # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)O"
        "PC": "{sn2}OC(CO{sn1})COP(=O)([O-])OCC[N+](C)(C)C", # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)([O-])OCC[N+](C)(C)C"
        "PE": "{sn2}OC(CO{sn1})COP(=O)(O)OCCN", # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCCN"
        "PI": "{sn2}OC(CO{sn1})COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O"
        "PS": "{sn2}OC(CO{sn1})COP(=O)(O)OCC(N)C(=O)O", # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCC(N)C(=O)O"
        "PG": "{sn2}OC(CO{sn1})COP(=O)(O)OCC(O)CO", # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCC(O)CO"
        "MG": "{sn1}OCC(O)CO", # f"{sn1_forward}OCC(O)CO"
        "DG": "{sn2}OC(CO)CO{sn1}", # f"{sn2_forward}OC(CO)CO{sn1_reverse}"
        "TG": "{sn2}OC(CO{sn3})CO{sn1}", # f"{sn2_forward}OC(CO{sn3_reverse})CO{sn1_reverse}"
        "Cer": "{sn1}{sn2}", # f"{sn1_forward}{sn2_reverse}" # fa chain specific parse
        "SM": "{sn1}C(O)C(COP(=O)([O-])OCC[N+](C)(C)C)N{sn2}", # f"{sn1_forward}C(O)C(COP(=O)([O-])OCC[N+](C)(C)C)N{sn2_reverse}"
        "BMP": "{sn2}OC(CO)COP(=O)(O)OCC(CO)O{sn1}", # f"{sn2_forward}OC(CO)COP(=O)(O)OCC(CO)O{sn1_reverse}"
        "CAR": "{sn1}OC(CC(=O)[O-])C[N+](C)(C)C", # f"{sn1_forward}OC(CC(=O)[O-])C[N+](C)(C)C"
        "LPA": "{sn1}OCC(O)COP(=O)(O)O", # f"{sn1_forward}OCC(O)COP(=O)(O)O"
        "LPC": "{sn1}OCC(O)COP(=O)([O-])OCC[N+](C)(C)C", # f"{sn1_forward}OCC(O)COP(=O)([O-])OCC[N+](C)(C)C"
        "LPE": "{sn1}OCC(O)COP(=O)(O)OCCN", # f"{sn1_forward}OCC(O)COP(=O)(O)OCCN"
        "LPI": "{sn1}OCC(O)COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", # f"{sn1_forward}OCC(O)COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O"
        "LPS": "{sn1}OCC(O)COP(=O)(O)OCC(N)C(=O)O", # f"{sn1_forward}OCC(O)COP(=O)(O)OCC(N)C(=O)O"
        "LPG": "{sn1}OCC(O)COP(=O)(O)OCC(O)CO", # f"{sn1_forward}OCC(O)COP(=O)(O)OCC(O)CO"
        "CerP": "{sn1}C(O)C(COP(=O)(O)O)N{sn2}", # f"{sn1_forward}C(O)C(COP(=O)(O)O)N{sn2_reverse}"
        "CE": "{sn1}OC1CCC2(C)C(=CCC3C2CCC2(C)C(C(C)CCCC(C)C)CCC32)C1", # f"{sn1_forward}OC1CCC2(C)C(=CCC3C2CCC2(C)C(C(C)CCCC(C)C)CCC32)C1"
        "CL": "{sn4}OCC(COP(=O)(O)OCC(O)COP(=O)(O)OCC(CO{sn1})O{sn2})O{sn3}", # f"{sn4_forward}OCC(COP(=O)(O)OCC(O)COP(=O)(O)OCC(CO{sn1_reverse})O{sn2_reverse})O{sn3_reverse}"
        # ...
    }
    FA_CHAIN_REVERSE = {
        "PA": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)O"
        "PC": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)([O-])OCC[N+](C)(C)C"
        "PE": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCCN"
        "PI": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O"
        "PS": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCC(N)C(=O)O"
        "PG": [True, False],    # f"{sn2_forward}OC(CO{sn1_reverse})COP(=O)(O)OCC(O)CO"
        "MG": [False],          # f"{sn1_forward}OCC(O)CO"
        "DG": [True, False],    # f"{sn2_forward}OC(CO)CO{sn1_reverse}"
        "TG": [True, False, True],  # f"{sn2_forward}OC(CO{sn3_reverse})CO{sn1_reverse}"
        "Cer": [False, True],   # f"{sn1_forward}{sn2_reverse}" # fa chain specific parse
        "SM": [False, True],    # f"{sn1_forward}C(O)C(COP(=O)([O-])OCC[N+](C)(C)C)N{sn2_reverse}"
        "BMP": [True, False],   # f"{sn2_forward}OC(CO)COP(=O)(O)OCC(CO)O{sn1_reverse}"
        "CAR": [False],         # f"{sn1_forward}OC(CC(=O)[O-])C[N+](C)(C)C"
        "LPA": [False],         # f"{sn1_forward}OCC(O)COP(=O)(O)O"
        "LPC": [False],         # f"{sn1_forward}OCC(O)COP(=O)([O-])OCC[N+](C)(C)C"
        "LPE": [False],         # f"{sn1_forward}OCC(O)COP(=O)(O)OCCN"
        "LPI": [False],         # f"{sn1_forward}OCC(O)COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O"
        "LPS": [False],         # f"{sn1_forward}OCC(O)COP(=O)(O)OCC(N)C(=O)O"
        "LPG": [False],         # f"{sn1_forward}OCC(O)COP(=O)(O)OCC(O)CO"
        "CerP": [False, True],  # f"{sn1_forward}C(O)C(COP(=O)(O)O)N{sn2_reverse}" 
        "CE": [False],          # f"{sn1_forward}OC1CCC2(C)C(=CCC3C2CCC2(C)C(C(C)CCCC(C)C)CCC32)C1"
        "CL": [True, True, True, False],    # f"{sn4_forward}OCC(COP(=O)(O)OCC(O)COP(=O)(O)OCC(CO{sn1_reverse})O{sn2_reverse})O{sn3_reverse}"
        # ...
    }    
    def __init__(
        self,
        repr_type: str = "lipidmaps",
        lipid_smiles_templates: Dict[str, str] = None, 
        fa_chain_reverse: Dict[str, bool] = None
    ):
        self.repr_type = repr_type
        if lipid_smiles_templates is not None:
            self.smiles_templates = lipid_smiles_templates
        else:
            self.smiles_templates = LipidSmilesGenerator.SMILES_TEMPLATES
        
        if fa_chain_reverse is not None:
            self.fa_chain_reverse = fa_chain_reverse
        else:
            self.fa_chain_reverse = LipidSmilesGenerator.FA_CHAIN_REVERSE
    
    def _parse_lipid_class(
        self,
        lipid_name: str,
        repr_type: str = "lipidmaps"
    ) -> str:
        """ Parse lipid class from lipid nomenclature.

        Parameters:
            lipid_name:
                Lipid nomenclature.
            repr_type:
                Type of lipid nomenclature.

        Returns:
            lipid_class:
                Lipid class of input lipid nomenclature.
        """
        match repr_type:
            case "lipidmaps":
                lipid_class = re.split(pattern=" ", string=lipid_name)[0]
            case "lipidsearch":
                # TODO
                pass
            case _:
                raise ValueError(
                    "lipid name representation error! \
                    Valid representation are `lipidmaps` and `lipidsearch`"
                )
    
        return lipid_class
    
    def _parse_fatty_acid(
        self,
        lipid_name: str,
        repr_type: str = "lipidmaps"
    ) -> str:
        """Parse fatty acids from lipid nomenclature.

        Parameters:
            lipid_name:
                Lipid nomenclature.
            repr_type:
                Type of lipid nomenclature.

        Returns:
            fatty_acid:
                List of fatty acids of input lipid nomenclature.
        """
        match repr_type:
            case "lipidmaps":
                fas_split = re.split(pattern=" ", string=lipid_name)
                if len(fas_split) > 1:
                    fas = fas_split[1]
                    fatty_acid_chains = re.split(pattern="[_/]", string=fas)
                else:
                    return None
            case "lipidsearch":
                # TODO
                pass
            case _:
                raise ValueError(
                    "lipid name representation error! \
                    Valid representation are `lipidmaps` and `lipidsearch`"
                )

        return fatty_acid_chains

    
    def _parse_fa_chain(
        self, 
        fa_chain: str
    ) -> Tuple[int, int, List[int]]:
        """Parse fatty acid chain notation."""
        new_fa_chain = re.sub(pattern=r";O\d+", repl="", string=fa_chain)
        new_fa_chain = re.sub(pattern=r"^(O-|P-)", repl="", string=new_fa_chain)
        new_fa_chain = re.sub(pattern=r"[ZE]", repl="", string=new_fa_chain)
        
        num_carbon_match = re.search(pattern=r"(\d+):", string=new_fa_chain)
        num_carbon = int(num_carbon_match.group(1)) if num_carbon_match else 0

        num_db_match = re.search(r":(\d+)", new_fa_chain)
        num_db = int(num_db_match.group(1)) if num_db_match else 0
        
        db_pos_match = re.search(r'\((\d+(?:,\d+)*)\)', new_fa_chain)
        db_positions = [int(p) for p in db_pos_match.group(1).split(",")] if db_pos_match else []

        return num_carbon, num_db, db_positions
    
    def _parse_hydroxyl(
        self, 
        fa_chain: str
    ) -> str:
        """Parse hydroxyl notation from fatty acid chain."""        
        hydroxyl_match = re.search(pattern=r";(O\d*)", string=fa_chain)
        if hydroxyl_match:
            hydroxyl = hydroxyl_match.group(1)
            return hydroxyl
        return ""    
    
    def _insert_after_n_carbon(
        self,
        smiles: str, 
        pos: int, 
        ins: str = "="
    ) -> str:
        """Insert a string after n-th carbon in the SMILES representation."""
        cnt = 0
        for idx, ch in enumerate(smiles):
            if ch == "C":
                cnt += 1
                if cnt == pos:
                    return smiles[:idx + 1] + ins + smiles[idx + 1:]
        return smiles
    
    def lipid_name_to_smiles(
        self, 
        lipid_name: str
    ) -> str:
        """"""
        lipid_class = self._parse_lipid_class(lipid_name=lipid_name, repr_type=self.repr_type)
        fa_chain_list = self._parse_fatty_acid(lipid_name=lipid_name, repr_type=self.repr_type)
        if fa_chain_list is None:
            return ""
        lipid_smiles = self.build_lipid_smiles(lipid_class, *fa_chain_list)
    
        return lipid_smiles
    
    def build_lipid_smiles(
        self, 
        lipid_class: str, 
        *fa_chains: str
    ) -> str:
        """
        exmple:
            fa_chain_list = ["16:0" , "18:1(9)"]
            smiles = build_lipid_smiles("PC", *fa_chain_list)
        
        """
        if lipid_class not in self.smiles_templates:
            raise KeyError(f"Unknown lipid class: {lipid_class}")
    
        reverses = self.fa_chain_reverse[lipid_class]
        if len(fa_chains) != len(reverses): 
            raise ValueError(f"{lipid_class} needs {len(reverses)} fatty acid chain, \
                         but given {len(fa_chains)}")
    
        parts = {
            f"sn{i+1}": self.fa_chain_to_smiles(
                            fa_chain=fa_chains[i], 
                            reverse=reverses[i], 
                            lipid_class=lipid_class
                        ) 
            for i in range(len(fa_chains))
        }
    
        return self.smiles_templates[lipid_class].format(**parts)
    
    def fa_chain_to_smiles(
        self,
        fa_chain: str, 
        reverse: bool = False,
        lipid_class: str = None,
    ) -> str:
        """Convert fatty acid chain notation to SMILES representation."""
        
        has_ozone = bool(re.search(r";O\d*", fa_chain))
    
        if has_ozone:
            if lipid_class in ["Cer"]:
                fa_smiles = self._ceramide_fa_chain_to_smiles(fa_chain=fa_chain)
            elif lipid_class in ["SM", "CerP"]:
                fa_smiles = self._other_sphingolipid_fa_chain_to_smiles(fa_chain=fa_chain)
        else: 
            if "P-" in fa_chain:
                fa_smiles = self._plasmanyl_ether_fa_chain_to_smiles(fa_chain=fa_chain, reverse=reverse)
            elif "O-" in fa_chain:
                fa_smiles = self._alkyl_ether_fa_chain_to_smiles(fa_chain=fa_chain, reverse=reverse)
            else:
                fa_smiles = self._normal_fa_chain_to_smiles(fa_chain=fa_chain, reverse=reverse)

        return fa_smiles

    def _normal_fa_chain_to_smiles(
        self, 
        fa_chain: str,
        reverse: bool = False
    ) -> str:
        """"""
        num_carbon, num_db, db_positions = self._parse_fa_chain(fa_chain=fa_chain)
        if num_db != len(db_positions):
            print("fa_chain:", fa_chain, "num_db:", num_db, "db_positions:", db_positions)
            raise ValueError("")
        init_smiles = "C" * (num_carbon - 1)
        
        if reverse:
            reverse_smiles = init_smiles
            for _, db_position in enumerate(db_positions):
                position = db_position - 1
                reverse_smiles = self._insert_after_n_carbon(
                                    smiles=reverse_smiles,
                                    pos=position,
                                    ins="="
                                )
            fa_smiles =  "C(=O)" + reverse_smiles
        else:
            forward_smiles = init_smiles
            for _, db_position in enumerate(db_positions):
                position = num_carbon - db_position
                forward_smiles = self._insert_after_n_carbon(
                                    smiles=forward_smiles,
                                    pos=position,
                                    ins="="
                                )
            fa_smiles = forward_smiles + "C(=O)"
        
        return fa_smiles

    def _alkyl_ether_fa_chain_to_smiles(
        self, 
        fa_chain: str,
        reverse: bool = False
    ) -> str:
        """"""
        num_carbon, num_db, db_positions = self._parse_fa_chain(fa_chain=fa_chain)
        if num_db != len(db_positions):
            raise ValueError("")
        init_smiles = "C" * (num_carbon - 1)
        
        if reverse:
            reverse_smiles = init_smiles
            for _, db_position in enumerate(db_positions):
                position = db_position - 1
                reverse_smiles = self._insert_after_n_carbon(
                                    smiles=reverse_smiles,
                                    pos=position,
                                    ins="="
                                )
            fa_smiles =  "C" + reverse_smiles
        else:
            forward_smiles = init_smiles
            for _, db_position in enumerate(db_positions):
                position = num_carbon - db_position
                forward_smiles = self._insert_after_n_carbon(
                                    smiles=forward_smiles,
                                    pos=position,
                                    ins="="
                                )
            fa_smiles = forward_smiles + "C"
        
        return fa_smiles

    def _plasmanyl_ether_fa_chain_to_smiles(
        self, 
        fa_chain: str,
        reverse: bool = False
    ) -> str:
        """"""
        num_carbon, num_db, db_positions = self._parse_fa_chain(fa_chain=fa_chain)
        if num_db != len(db_positions):
            raise ValueError("")
        init_smiles = "C" * (num_carbon - 2)
        
        if reverse:
            reverse_smiles = init_smiles
            for _, db_position in enumerate(db_positions):
                position = db_position - 2
                reverse_smiles = self._insert_after_n_carbon(
                                    smiles=reverse_smiles,
                                    pos=position,
                                    ins="="
                                )
            fa_smiles =  "C=C" + reverse_smiles
        else:
            forward_smiles = init_smiles
            for _, db_position in enumerate(db_positions):
                position = num_carbon - db_position
                forward_smiles = self._insert_after_n_carbon(
                                    smiles=forward_smiles,
                                    pos=position,
                                    ins="="
                                )
            fa_smiles = forward_smiles + "C=C"

        return fa_smiles
    
    def _ceramide_fa_chain_to_smiles(
        self, 
        fa_chain: str
    ) -> str:
        """"""    
        num_carbon, num_db, db_positions = self._parse_fa_chain(fa_chain=fa_chain)
        hydroxyl = self._parse_hydroxyl(fa_chain=fa_chain)
        if num_db != len(db_positions):
            raise ValueError("")
        init_smiles = "C" * (num_carbon - 3)
        
        match hydroxyl:
            case "O":
                backbone_smiles = "C(O)C(C)N"
                forward_smiles = init_smiles
            case "O2":
                backbone_smiles = "C(O)C(CO)N"
                forward_smiles = init_smiles
            case "O3":
                backbone_smiles = "C(O)C(CO)N" # C(O)C(O)C(CO)N
                forward_smiles = init_smiles + "(O)"
            
        for _, db_position in enumerate(db_positions):
            position = num_carbon - db_position
            forward_smiles = self._insert_after_n_carbon(
                                    smiles=forward_smiles,
                                    pos=position,
                                    ins="="
                                )
        fa_smiles = forward_smiles + backbone_smiles
    
        return fa_smiles

    def _other_sphingolipid_fa_chain_to_smiles(
        self, 
        fa_chain: str
    ) -> str:
        """"""
        num_carbon, num_db, db_positions = self._parse_fa_chain(fa_chain=fa_chain)
        hydroxyl = self._parse_hydroxyl(fa_chain=fa_chain)
        if num_db != len(db_positions):
            raise ValueError("")
        init_smiles = "C" * (num_carbon - 3)
        
        match hydroxyl:
            case "O":
                return "" 
            case "O2":
                forward_smiles = init_smiles
            case "O3":
                forward_smiles = init_smiles + "(O)"
        
        for _, db_position in enumerate(db_positions):
            position = num_carbon - db_position
            forward_smiles = self._insert_after_n_carbon(
                                    smiles=forward_smiles,
                                    pos=position,
                                    ins="="
                                )
        fa_smiles = forward_smiles
        
        return fa_smiles
        
    
if __name__ == "__main__":
    in_file = "clean_data/lipidmaps_chain_shorthan.csv"
    out_file = "clean_data/lipidmaps_chain_shorthand_with_smiles.csv"
    
    df_lipid = pd.read_csv(in_file)
    smiles_generator = LipidSmilesGenerator(repr_type="lipidmaps")
    
    df_lipid["smiles"] = df_lipid["db_shorthand"].progress_apply(
        lambda lipid_name: smiles_generator.lipid_name_to_smiles(lipid_name=lipid_name)
    )
    df_lipid = df_lipid[df_lipid["smiles"].astype(bool)]
    df_lipid.to_csv(out_file, index=False)
    
    
"""Mass spectrometry data and lipid database loading utilities.

This module provides classes to load and parse:
1. Lipid ozonolysis product reference databases (JSON format)
2. MS1 peak tables from MZmine or MS-DIAL
3. MS2 mass spectra from MGF or MSP files

It standardizes data formats for downstream lipid identification and analysis.

Example:
    >>> db_loader = DatabaseLoader("lipid_database.json")
    >>> ref_df = db_loader.load_reference_database()
    >>> data_loader = DataLoader("ms1.csv", "ms2.mgf")
    >>> ms1_df = data_loader.load_ms1_peak_table()
    >>> ms2_spectra = data_loader.load_ms2_spectrum()

Classes:
    ReferenceLipid: Data model for lipid reference information.
    DatabaseLoader: Loader for ozonolysis product reference databases.
    DataLoader: Unified loader for MS1 peak tables and MS2 spectral data.
"""

from typing import List
import json

import pandas as pd
from matchms import Spectrum
from matchms.importing import load_from_mgf, load_from_msp


class ReferenceLipid:
    """Data model class to store lipid reference information.

    Attributes:
        primary_annotated_name: Primary annotated name of the lipid.
        structure_annotated_name: Detailed structure annotated name.
        fragment_mz_delta: List of m/z delta (mass shift) values for fragments.
        fragment_fatty_acid_pos: List of delta positions of corresponding fatty acids.
        fragment_types: List of fragment types (e.g., 'aldehyde' or 'criegee').
    """
    
    def __init__(
        self,
        primary_annotated_name: str,
        structure_annotated_name: str,
        fragment_mz_delta: List[float],
        fragment_fatty_acid_pos: List[int],
        fragment_types: List[str]
    ):
        """Initializes ReferenceLipid instance with fragment and naming data.

        Args:
            primary_annotated_name: Primary name for lipid annotation.
            structure_annotated_name: Structural shorthand name.
            fragment_mz_delta: List of fragment mass delta values.
            fragment_fatty_acid_pos: List of fatty acid delta positions.
            fragment_types: List of fragment classification types.
        """
        self.primary_annotated_name = primary_annotated_name
        self.structure_annotated_name = structure_annotated_name
        self.fragment_mz_delta = fragment_mz_delta
        self.fragment_fatty_acid_pos = fragment_fatty_acid_pos
        self.fragment_types = fragment_types


class DatabaseLoader:
    """Loader class for reading and parsing ozonolysis product reference database.

    Attributes:
        ozid_database_file: Path to the JSON ozonolysis product database file.
    """
    
    def __init__(
        self,
        ozid_database_file: str
    ):
        """Initializes database loader with target file path.

        Args:
            ozid_database_file: Path to the JSON reference database file.
        """
        self.ozid_database_file = ozid_database_file

    def load_reference_database(self) -> pd.DataFrame:
        """Loads and returns formatted lipid reference database.

        Returns:
            pd.DataFrame: DataFrame containing standardized lipid reference data.
        """
        df_reference = self._load_ozonolysis_prodcut_database()
        return df_reference

    def _load_ozonolysis_prodcut_database(self) -> pd.DataFrame:

        """Loads and processes ozonolysis product data from JSON file.

        Returns:
            pd.DataFrame: DataFrame with renamed and structured database entries.
        """
        with open(self.ozid_database_file, "r") as f:
            delta_mass_list = json.load(f)
        df = pd.DataFrame().from_records(delta_mass_list)
        df = df.rename(columns={
                "chain_shorthand": "primary_annotated_name",
                "db_shorthand": "structure_annotated_name",
                "smiles": "smiles",
                "delta_mass": "fragment_mz_delta",
                "fatty_acid_pos": "fragment_fatty_acid_pos",
                "fragment_type": "fragment_types"
            })
        return df


class DataLoader:
    """Unified loader for MS1 peak tables and MS2 spectral data.

    Supports multiple input formats:
    MS1: MZmine, MS-DIAL
    MS2: MGF, MSP

    Attributes:
        ms1_file: Path to MS1 peak table file.
        ms2_file: Path to MS2 spectrum file.
        ms1_file_type: Type of MS1 input (default: 'mzmine').
        ms2_file_type: Type of MS2 input (default: 'mgf').
    """
    
    def __init__(
        self,
        ms1_file: str,
        ms2_file: str,
        ms1_file_type: str = "mzmine",
        ms2_file_type: str = "mgf"
    ):
        """Initializes mass spec data loader with input files and formats.

        Args:
            ms1_file: Path to MS1 peak table file.
            ms2_file: Path to MS2 spectrum file.
            ms1_file_type: Format of MS1 file ('mzmine' or 'msdial').
            ms2_file_type: Format of MS2 file ('mgf' or 'msp').
        """
        self.ms1_file = ms1_file
        self.ms2_file = ms2_file
        self.ms1_file_type = ms1_file_type
        self.ms2_file_type = ms2_file_type

    def load_ms1_peak_table(self) -> pd.DataFrame:
        """Loads and parses MS1 peak table based on configured input type.

        Returns:
            pd.DataFrame: Standardized DataFrame of MS1 features and intensities.

        Raises:
            ValueError: If unsupported ms1_file_type is provided.
        """
        match self.ms1_file_type:
            case "mzmine":
                df_ms1 = self._load_mzmine_peak_table()
            case "msdial":
                df_ms1 = self._load_msdial_peak_table()
            case _:
                raise ValueError(f"Unsupported ms1_file_type: {self.ms1_file_type}")

        return df_ms1

    def _load_mzmine_peak_table(self) -> pd.DataFrame:
        """Loads and cleans peak table output from MZmine software.

        Returns:
            pd.DataFrame: Formatted MS1 table with feature IDs, m/z, RT, and intensities.
        """
        df = pd.read_csv(self.ms1_file, low_memory=False)
        feature_cols = [
            "id", "rt", "mz",
            "spectral_db_matches:precursor_mz",
            "spectral_db_matches:mol_formula",
            "spectral_db_matches:ion_adduct", 
            "spectral_db_matches:compound_name", 
        ]
        data_cols = df.columns[df.columns.str.match("datafile.*area")].to_list()
        df_data = df[data_cols].fillna(0.0)

        df_vars = df[feature_cols]
        df_vars = df_vars.sort_values(["rt", "mz"], ascending=[True, True])
        df_vars.columns = df_vars.columns.str.replace(r"spectral_db_matches:", "")
        df_vars = df_vars.rename(columns={
                                    "precursor_mz": "mz_calculated",
                                })
        df_vars = df_vars.rename(columns={
                                    "id": "feature_id",
                                    "mz": "precursor_mz",
                                    "ion_adduct": "adduct",
                                    "compound_name": "primary_annotated_name"
                                })
        df_vars = df_vars.fillna("")
        
        ms1_table = pd.concat([df_vars, df_data],axis=1)
        ms1_table["feature_id"] = ms1_table["feature_id"].astype(str)
        
        return ms1_table

    def _load_msdial_peak_table(self) -> pd.DataFrame:
        """Placeholder for loading MS-DIAL output peak table (not implemented).

        Returns:
            pd.DataFrame: Formatted MS1 table from MS-DIAL.
        """
        # TODO
        pass

    def load_ms2_spectrum(self) -> List[Spectrum]:
        """Loads MS2 spectra from file based on configured format.

        Returns:
            List[Spectrum]: List of matchms Spectrum objects.

        Raises:
            ValueError: If unsupported ms2_file_type is provided.
        """
        match self.ms2_file_type:
            case "mgf":
                ms2_spectra = self._load_spectrum_from_mgf()
            case "msp":
                ms2_spectra = self._load_spectrum_from_msp()
            case _:
                raise ValueError(f"Unsupported ms2_file_type: {self.ms2_file_type}")

        return ms2_spectra

    def _load_spectrum_from_mgf(self) -> List[Spectrum]:
        """Loads MS2 spectra from MGF format file using matchms.

        Returns:
            List[Spectrum]: List of parsed Spectrum objects.
        """
        spectra = list(load_from_mgf(self.ms2_file))
        return spectra

    def _load_spectrum_from_msp(self) -> List[Spectrum]:
        """Loads MS2 spectra from MSP format file using matchms.

        Returns:
            List[Spectrum]: List of parsed Spectrum objects.
        """
        spectra = list(load_from_msp(self.ms2_file))
        return spectra




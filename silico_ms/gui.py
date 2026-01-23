from typing import List, Dict
import os


import pandas as pd
from matchms import Spectrum
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt6.QtWidgets import (
    QMainWindow, QWidget,
    QLineEdit, QLabel, QComboBox,
    QDialog, QFileDialog,
    QPushButton,
    QVBoxLayout, QHBoxLayout,
    QMessageBox, QCheckBox,
    QTabWidget, 
    QTableWidget, QTableWidgetItem
)
from PyQt6.QtGui import (
    QAction,
    QIntValidator, QDoubleValidator
)
from PyQt6.QtCore import (
    QObject, pyqtSignal
)


from silico_ms.data_utils import (
    ReferenceLipid, DataLoader, DatabaseLoader,
    STRUCTURE_FILE_PATH, OZID_FILE_PATH
)
from silico_ms.spectrum_utils import SpectrumFeatureProcesser, split_features
from silico_ms.plot_utils import spectrum_plot, none_spectrum_plot, fragment_network_plot
from silico_ms.algorithm import FragmentNetowrk, LipidAnnotator



class DataImportUI(QObject):
    
    project_name = pyqtSignal(str)
    ms1_peak_table = pyqtSignal(pd.DataFrame)
    ms2_spectra = pyqtSignal(object) # List[Spectrum]
    reference_database = pyqtSignal(object) # List[ReferenceLipid]

    def __init__(
        self,
        parent: QObject = None
    ) -> None:
        super().__init__(parent)

    def create_menu_actions(
        self,
        parent: QObject = None
    ) -> None:
        self.import_data_act = QAction("Import raw data", parent)
        self.import_database_act = QAction("Import reference database", parent)

        self.import_data_act.triggered.connect(self.import_raw_data)
        self.import_database_act.triggered.connect(self.import_reference_database)

    def import_raw_data(self) -> None:
        raw_data_import = RawDataImportDialog(self.parent())

        raw_data_import.ms1_peak_table.connect(self.on_ms1_peak_table_loaded)
        raw_data_import.ms2_spectra.connect(self.on_ms2_spectra_loaded)
        raw_data_import.project_name.connect(self.on_project_name_loaded)

        raw_data_import.show()

    def import_reference_database(self):
        database_import = ReferenceDatabaseImportDialog(self.parent())

        database_import.reference_database.connect(self.on_reference_database_loaded)

        database_import.show()

    def on_ms1_peak_table_loaded(
        self,
        peak_table: pd.DataFrame
    ):
        self.ms1_peak_table.emit(peak_table)

    def on_ms2_spectra_loaded(
        self,
        spectra: List[Spectrum]
    ):
        self.ms2_spectra.emit(spectra)

    def on_reference_database_loaded(
        self,
        database: List[ReferenceLipid]
    ):
        self.reference_database.emit(database)
        
    def on_project_name_loaded(
        self,
        name: str
    ):
        self.project_name.emit(name)


class RawDataImportDialog(QDialog):

    project_name = pyqtSignal(str)
    ms1_peak_table = pyqtSignal(pd.DataFrame)
    ms2_spectra = pyqtSignal(object) # List[Spectrum]

    def __init__(
        self,
        parent: QObject = None
    ) -> None:
        super().__init__(parent)

        self.init_ui()

    def init_ui(self):
        self.ms1_file_line_edit = QLineEdit()
        self.ms2_file_line_edit = QLineEdit()
        self.project_name_line_edit = QLineEdit()
        self.project_name_line_edit.setPlaceholderText("project001")
        self.ms1_file_line_edit.setReadOnly(True)
        self.ms2_file_line_edit.setReadOnly(True)

        self.ms1_file_type_combo_box = QComboBox()
        self.ms1_file_type_combo_box.addItem("Mzmine", "mzmine")
        self.ms1_file_type_combo_box.addItem("MS-DIAL", "msdial")
        self.ms2_file_type_combo_box = QComboBox()
        self.ms2_file_type_combo_box.addItem(".msp", "msp")
        self.ms2_file_type_combo_box.addItem(".mgf", "mgf")

        self.ms1_file_browse_button = QPushButton("Select MS1 peak table file")
        self.ms2_file_browse_button = QPushButton("Select MS2 spectra file")
        self.ms1_file_browse_button.clicked.connect(self.import_ms1_file)
        self.ms2_file_browse_button.clicked.connect(self.import_ms2_file)

        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")
        self.ok_button.clicked.connect(self.load_all_data)
        self.cancel_button.clicked.connect(self.reject)

        layout = QVBoxLayout()

        ms1_file_layout = QHBoxLayout()
        ms1_file_layout.addWidget(QLabel("MS1 peak table file name"))
        ms1_file_layout.addWidget(self.ms1_file_line_edit)
        ms1_file_layout.addWidget(self.ms1_file_type_combo_box)
        ms1_file_layout.addWidget(self.ms1_file_browse_button)

        ms2_file_layout = QHBoxLayout()
        ms2_file_layout.addWidget(QLabel("MS2 spectra file name"))
        ms2_file_layout.addWidget(self.ms2_file_line_edit)
        ms2_file_layout.addWidget(self.ms2_file_type_combo_box)
        ms2_file_layout.addWidget(self.ms2_file_browse_button)

        project_layout = QHBoxLayout()
        project_layout.addWidget(QLabel("Project name"))
        project_layout.addWidget(self.project_name_line_edit)

        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(self.ok_button)
        buttons_layout.addWidget(self.cancel_button)

        layout.addLayout(ms1_file_layout)
        layout.addLayout(ms2_file_layout)
        layout.addLayout(project_layout)
        layout.addLayout(buttons_layout)

        self.setLayout(layout)
        self.setWindowTitle("Import raw data")
        self.resize(600, 400)

    def import_ms1_file(self):
        file_path = self._select_file()

        if file_path:
            self.ms1_file_line_edit.setText(file_path)
            self.ms1_file_path = file_path

    def import_ms2_file(self):
        file_path = self._select_file()

        if file_path:
            self.ms2_file_line_edit.setText(file_path)
            self.ms2_file_path = file_path

    def _select_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
                            parent=self.parent(),
                            caption="Select files",
                            directory="",
                            filter="All files (*);;Csv file (*.csv)"
                        )

        return file_path if file_path else None

    def get_base_name_from_path(self, path: str) -> str:
        """"""
        file_name = os.path.basename(path)
        return os.path.splitext(file_name)[0]

    def load_all_data(self):
        """"""
        ms1_file_type = self.ms1_file_type_combo_box.currentData()
        ms2_file_type = self.ms2_file_type_combo_box.currentData()
        task_name = self.get_base_name_from_path(path=self.ms1_file_path)

        data_loader = DataLoader(
                        ms1_file=self.ms1_file_path,
                        ms2_file=self.ms2_file_path,
                        ms1_file_type=ms1_file_type,
                        ms2_file_type=ms2_file_type
                    )

        peak_table = data_loader.load_ms1_peak_table()
        spectra = data_loader.load_ms2_spectrum()

        self.ms1_peak_table.emit(peak_table)
        self.ms2_spectra.emit(spectra)
        self.project_name.emit(task_name)

        self.accept()
        QMessageBox.information(
            self,
            "Import data",
            "Loading data Done!",
            buttons=QMessageBox.StandardButton.Ok
        )


class ReferenceDatabaseImportDialog(QDialog):

    reference_database = pyqtSignal(object) # List[ReferenceLipid]

    def __init__(
        self,
        parent: QObject = None
    ) -> None:
        super().__init__(parent)

        self.init_ui()

    def init_ui(self):
        """"""
        self.ozid_file_line_edit = QLineEdit()
        self.structure_file_line_edit = QLineEdit()
        self.ozid_file_line_edit.setReadOnly(True)
        self.structure_file_line_edit.setReadOnly(True)

        self.use_customize_check_box = QCheckBox()
        self.use_customize_check_box.setChecked(False)
        self.use_customize_check_box.stateChanged.connect(self.on_checkbox_changed)
        self.ozid_file_browse_button = QPushButton("Select ozonolysis prodcut database file")
        self.structure_file_browse_button = QPushButton("Select lipid structure database file")
        self.ozid_file_browse_button.clicked.connect(self.import_ozid_file)
        self.structure_file_browse_button.clicked.connect(self.import_structure_file)

        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")
        self.ok_button.clicked.connect(self.load_all_database)
        self.cancel_button.clicked.connect(self.reject)

        main_layout = QVBoxLayout()

        use_customize_layout = QHBoxLayout()
        use_customize_layout.addWidget(QLabel("Use customize"))
        use_customize_layout.addWidget(self.use_customize_check_box)

        ozid_file_layout = QHBoxLayout()
        ozid_file_layout.addWidget(QLabel("Onolysis prodcut database file name"))
        ozid_file_layout.addWidget(self.ozid_file_line_edit)
        ozid_file_layout.addWidget(self.ozid_file_browse_button)

        structure_file_layout = QHBoxLayout()
        structure_file_layout.addWidget(QLabel("Structure database file name"))
        structure_file_layout.addWidget(self.structure_file_line_edit)
        structure_file_layout.addWidget(self.structure_file_browse_button)

        self.customize_widget = QWidget()
        customize_layout = QVBoxLayout()
        customize_layout.addLayout(ozid_file_layout)
        customize_layout.addLayout(structure_file_layout)
        self.customize_widget.setLayout(customize_layout)
        self.customize_widget.setVisible(False)

        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(self.ok_button)
        buttons_layout.addWidget(self.cancel_button)

        main_layout.addLayout(use_customize_layout)
        main_layout.addWidget(self.customize_widget)
        main_layout.addLayout(buttons_layout)

        self.setLayout(main_layout)
        self.setWindowTitle("Import reference database")
        self.resize(600, 400)

    def on_checkbox_changed(self):
        if self.use_customize_check_box.isChecked():
            self.customize_widget.setVisible(True)
        else:
            self.customize_widget.setVisible(False)

    def import_ozid_file(self):
        file_path = self._select_file()

        if file_path:
            self.ozid_file_line_edit.setText(file_path)
            self.ozid_file_path = file_path

    def import_structure_file(self):
        file_path = self._select_file()

        if file_path:
            self.structure_file_line_edit.setText(file_path)
            self.structure_file_path = file_path

    def _select_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
                            parent=self.parent(),
                            caption="Select files",
                            directory="",
                            filter="All files (*);;Csv file (*.csv)"
                        )

        return file_path if file_path else None

    def load_all_database(self):
        if self.use_customize_check_box.isChecked():
            if  not (hasattr(self, "ozid_file_path") and
                    hasattr(self, "structure_file_path")):
                raise ValueError("Error, No `ozid_file_path` or `structure_file_path`")

            database_loader = DatabaseLoader(
                            ozid_database_file=self.ozid_file_path,
                            structure_database_file=self.structure_file_path
                        )
        else:
            database_loader = DatabaseLoader(
                            ozid_database_file=OZID_FILE_PATH,
                            structure_database_file=STRUCTURE_FILE_PATH
                        )

        database = database_loader.load_reference_database()
        self.reference_database.emit(database)

        self.accept()
        QMessageBox.information(
            self,
            "Import reference database",
            "Loading reference database done!",
            buttons=QMessageBox.StandardButton.Ok
        )


class LipidAnnotationUI(QObject):

    peak_table_annotated = pyqtSignal(pd.DataFrame)
    fragment_newtork_dict = pyqtSignal(dict) # Dict[str, FragmentNetwork]

    def __init__(self, parent: QObject = None):
        super().__init__(parent)

    def create_menu_actions(
        self,
        parent: QObject = None
    ) -> None:
        self.double_bond_annotation_act = QAction("C=C double bond annotation", parent)

        self.double_bond_annotation_act.triggered.connect(self.run_double_bond_annotation)

    def run_double_bond_annotation(self):
        double_bond_annotater = DoulbeBondAnnoatatorialog(self.parent())

        double_bond_annotater.peak_table_annotated.connect(self.on_double_bond_annotated)
        double_bond_annotater.fragment_newtork_dict.connect(self.on_fragment_network_generated)

        double_bond_annotater.show()

    def on_double_bond_annotated(
        self,
        peak_table: pd.DataFrame
    ):
        self.peak_table_annotated.emit(peak_table)

    def on_fragment_network_generated(
        self,
        newtowrk_dict: Dict[str, FragmentNetowrk]
    ):
        self.fragment_newtork_dict.emit(newtowrk_dict)


class DoulbeBondAnnoatatorialog(QDialog):

    peak_table_annotated = pyqtSignal(pd.DataFrame)
    fragment_newtork_dict = pyqtSignal(dict) # Dict[str, FragmentNetwork]

    def __init__(
        self,
        parent: QObject = None
    ) -> None:
        super().__init__(parent)

        if parent:
            parent.ms1_peak_table.connect(self.on_ms1_peak_table_loaded)
            parent.ms2_spectra.connect(self.on_ms2_spectra_loaded)

        self.init_ui()

    def init_ui(self) -> None:
        self.rt_tol_line_edit = QLineEdit()
        self.rt_tol_line_edit.setValidator(QDoubleValidator())
        self.rt_tol_mode_combo_box = QComboBox()
        self.rt_tol_mode_combo_box.addItem("absolute (min)", "absolute")
        self.rt_tol_mode_combo_box.addItem("relative", "relative")
        self.mz_tol_line_edit = QLineEdit()
        self.mz_tol_line_edit.setValidator(QDoubleValidator())
        self.mz_tol_mode_combo_box = QComboBox()
        self.mz_tol_mode_combo_box.addItem("Da", "Da")
        self.mz_tol_mode_combo_box.addItem("ppm", "ppm")

        self.rt_weight_line_edit = QLineEdit()
        self.rt_weight_line_edit.setValidator(QDoubleValidator())
        self.precurosr_mz_weight_line_edit = QLineEdit()
        self.precurosr_mz_weight_line_edit.setValidator(QDoubleValidator())
        self.spec_weight_line_edit = QLineEdit()
        self.spec_weight_line_edit.setValidator(QDoubleValidator())

        self.spec_similarity_type_combo_box = QComboBox()
        self.spec_similarity_type_combo_box.addItem("NIST-LC", "NIST-LC")
        self.spec_similarity_type_combo_box.addItem("NIST-GC", "NIST-GC")
        self.spec_similarity_type_combo_box.addItem("SQRT", "SQRT")
        self.spec_similarity_type_combo_box.addItem("MassBank", "MassBank")
        self.spec_similarity_type_combo_box.addItem("ModifiedCosine", "ModifiedCosine")
        self.spec_similarity_type_combo_box.addItem("NeutralLossesCosine", "NeutralLossesCosine")
        self.spec_similarity_type_combo_box.addItem("CosineHungarian", "CosineHungarian")
        self.spec_similarity_type_combo_box.addItem("None", "None")
        self.score_threshold_line_edit = QLineEdit()
        self.score_threshold_line_edit.setValidator(QDoubleValidator())
        self.top_n_line_edit = QLineEdit()
        self.top_n_line_edit.setValidator(QIntValidator())

        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")
        self.ok_button.clicked.connect(self.run_double_bond_annotation)
        self.cancel_button.clicked.connect(self.reject)

        main_layout = QVBoxLayout()

        rt_layout = QHBoxLayout()
        rt_layout.addWidget(QLabel("RT tolerance"))
        rt_layout.addWidget(self.rt_tol_line_edit)
        rt_layout.addWidget(QLabel("RT tolerance mode"))
        rt_layout.addWidget(self.rt_tol_mode_combo_box)

        mz_layout = QHBoxLayout()
        mz_layout.addWidget(QLabel("m/z tolerance"))
        mz_layout.addWidget(self.mz_tol_line_edit)
        mz_layout.addWidget(QLabel("m/z tolerance mode"))
        mz_layout.addWidget(self.mz_tol_mode_combo_box)

        weight_layout = QHBoxLayout()
        weight_layout.addWidget(QLabel("RT weight: "))
        weight_layout.addWidget(self.rt_weight_line_edit)
        weight_layout.addWidget(QLabel("m/z weight: "))
        weight_layout.addWidget(self.precurosr_mz_weight_line_edit)
        weight_layout.addWidget(QLabel("MS2 spectra weight: "))
        weight_layout.addWidget(self.spec_weight_line_edit)

        other_layout = QHBoxLayout()
        other_layout.addWidget(QLabel("MS2 spectrum similarity type: "))
        other_layout.addWidget(self.spec_similarity_type_combo_box)
        other_layout.addWidget(QLabel("Score threshold: "))
        other_layout.addWidget(self.score_threshold_line_edit)
        other_layout.addWidget(QLabel("Top N isomer: "))
        other_layout.addWidget(self.top_n_line_edit)

        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(self.ok_button)
        buttons_layout.addWidget(self.cancel_button)

        main_layout.addLayout(rt_layout)
        main_layout.addLayout(mz_layout)
        main_layout.addLayout(weight_layout)
        main_layout.addLayout(other_layout)
        main_layout.addLayout(buttons_layout)

        self.setLayout(main_layout)
        self.setWindowTitle("Lipid C=C double bond annotation")
        self.resize(600, 400)

    def run_double_bond_annotation(self):
        """"""
        ms1_peak_table = self.get_ms1_peak_table()
        ms2_spectra = self.get_ms2_spectra()
        rt_tol = float(self.rt_tol_line_edit.text())
        rt_tol_mode = self.rt_tol_mode_combo_box.currentData()
        mz_tol = float(self.mz_tol_line_edit.text())
        mz_tol_mode = self.mz_tol_mode_combo_box.currentData()
        rt_weight = float(self.rt_weight_line_edit.text())
        precusro_mz_weight = float(self.precurosr_mz_weight_line_edit.text())
        spec_weight = float(self.spec_weight_line_edit.text())
        ms2_spectrum_similarity_type = self.spec_similarity_type_combo_box.currentData()
        reference_database = self.get_reference_database()
        score_threshold = float(self.score_threshold_line_edit.text())
        top_n = int(self.top_n_line_edit.text())

        candidate_features, auxiliary_features = split_features(
                                                    ms1_peak_table=ms1_peak_table,
                                                    ms2_spectra=ms2_spectra
                                                )
        feature_processer  = SpectrumFeatureProcesser(
                                rt_tol=rt_tol,
                                rt_tol_mode=rt_tol_mode,
                                mz_tol=mz_tol,
                                mz_tol_mode=mz_tol_mode,
                                rt_weight=rt_weight,
                                precusro_mz_weight=precusro_mz_weight,
                                spec_weight=spec_weight,
                                ms2_spectrum_similarity_type=ms2_spectrum_similarity_type
                            )
        lipid_annotator = LipidAnnotator(
                            reference_database=reference_database,
                            feature_processer=feature_processer,
                            score_threshold=score_threshold,
                            top_n=top_n
                        )

        results = lipid_annotator.annotate_features(
                        candidate_features=candidate_features,
                        auxiliary_features=auxiliary_features
                    )
        network_dict = lipid_annotator.get_network_dict()

        self.df_result.emit(results)
        self.fragment_newtork_dict.emit(network_dict)

        self.accept()
        QMessageBox.information(
            self,
            "Annotate C=C double bond",
            "Annotate C=C double bond done!",
            buttons=QMessageBox.StandardButton.Ok
        )

    def get_ms1_peak_table(self) -> pd.DataFrame:
        if not hasattr(self, "ms1_peak_table"):
            raise ValueError("No MS1 peak table")

        return self.ms1_peak_table

    def get_ms2_spectra(self) -> List[Spectrum]:
        if not hasattr(self, "ms2_spectra"):
            raise ValueError("No MS2 spectra")

        return self.ms2_spectra

    def get_reference_database(self) -> List[ReferenceLipid]:
        if not hasattr(self, "reference_database"):
            raise ValueError("No reference database")

        return self.reference_database

    def on_ms1_peak_table_loaded(
        self,
        ms1_peak_table: pd.DataFrame
    ):
        self.ms1_peak_table = ms1_peak_table

    def on_ms2_spectra_loaded(
        self,
        ms2_spectra: List[Spectrum]
    ):
        self.ms2_spectra = ms2_spectra

    def on_reference_database_loaded(
        self,
        reference_database: List[ReferenceLipid]
    ):
        self.reference_database = reference_database


class DataExportUI(QObject):
    """
    """
    def __init__(
        self,
        parent: QObject = None
    ):
        super().__init__(parent)

    def create_menu_actions(
        self,
        parent: QObject = None
    ) -> None:
        self.data_export_act = QAction("Data export", parent)

        self.data_export_act.triggered.connect(self.export_data)

    def export_data(self):
        data_export = DataExportDialog(self.parent())

        data_export.show()


class DataExportDialog(QDialog):
    """
    """
    def __init__(
        self,
        parent: QObject = None
    ):
        super().__init__(parent)

        if parent:
            parent.peak_table_annotated.connect(self.on_results_loaded)

        self.init_ui()

    def init_ui(self):

        self.output_file_line_edit = QLineEdit()
        self.output_file_line_edit.setReadOnly(True)
        self.output_foler_browse_button = QPushButton("Select Output folder")
        self.output_foler_browse_button.clicked.connect(self.select_output_folder)

        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")
        self.ok_button.clicked.connect(self.export_results)
        self.cancel_button.clicked.connect(self.reject)

        main_layout = QVBoxLayout()

        output_file_layout = QHBoxLayout()
        output_file_layout.addWidget(QLabel("RT tolerance"))
        output_file_layout.addWidget(self.output_file_line_edit)
        output_file_layout.addWidget(self.output_foler_browse_button)

        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(self.ok_button)
        buttons_layout.addWidget(self.cancel_button)

        main_layout.addLayout(output_file_layout)
        main_layout.addLayout(buttons_layout)

        self.setLayout(main_layout)
        self.setWindowTitle("Results export")
        self.resize(600, 400)

    def export_results(self):
        if not hasattr(self, "folder_path"):
            raise ValueError("No output foler")
        if not hasattr(self, "df_output"):
            raise ValueError("No output dataframe")

        out_file = os.path.join(self.folder_path, "annoatated_feature.csv")
        self.df_output.to_csv(out_file, index=None)

    def select_output_folder(self):
        folder_path = QFileDialog.getExistingDirectory(
            self,
            "Select folder",
            "",
            QFileDialog.Option.ShowDirsOnly
        )
        if folder_path:
            self.folder_path = folder_path

    def on_results_loaded(
        self,
        df_results: pd.DataFrame
    ):
        self.df_output = df_results


class VisualizationUI(QObject):

    def __init__(
        self,
        parent: QObject = None
    ):
        super().__init__(parent)

    def create_menu_actions(
        self,
        parent: QObject = None
    ) -> None:
        pass


class MainWindow(QMainWindow):

    project_name = pyqtSignal(str)
    ms1_peak_table = pyqtSignal(pd.DataFrame)
    ms2_spectra = pyqtSignal(object) # List[Spectrum]
    reference_database = pyqtSignal(object) # List[ReferenceLipid]
    peak_table_annotated = pyqtSignal(pd.DataFrame)
    fragment_newtork_dict = pyqtSignal(dict) # Dict[str, FragmentNetwork]

    def __init__(self):
        super().__init__()

        self.data_import = DataImportUI(self)
        self.lipid_annotation= LipidAnnotationUI(self)
        self.data_export = DataExportUI(self)
        self.visualization = VisualizationUI(self)

        self.init_ui()

    def init_ui(self):
        """"""
        self.init_menu()
        self.init_widget()
        
        self.setGeometry(300, 300, 1600, 1200)
        self.setWindowTitle("Lintomics")
        self.show()

    def init_menu(self):
        self.menu_bar = self.menuBar()
        self.data_import_menu = self.menu_bar.addMenu("Data import")
        self.lipid_annotation_menu = self.menu_bar.addMenu("Lipid annotation")
        self.data_export_menu = self.menu_bar.addMenu("Results export")
        self.visualization_menu = self.menu_bar.addMenu("Visualization")

        self.data_import.create_menu_actions(self)
        self.data_import_menu.addAction(self.data_import.import_data_act)
        self.data_import_menu.addAction(self.data_import.import_database_act)
        self.data_import.ms1_peak_table.connect(self.on_raw_ms1_peak_table_loaded)
        self.data_import.ms2_spectra.connect(self.on_ms2_spectra_loaded)
        self.data_import.reference_database.connect(self.on_reference_database_loaded)
        self.data_import.project_name.connect(self.on_project_name_loaded)
        
        self.lipid_annotation.create_menu_actions(self)
        self.lipid_annotation_menu.addAction(self.lipid_annotation.double_bond_annotation_act)
        self.lipid_annotation.peak_table_annotated.connect(self.on_annotated_peak_table_loaded)
        self.lipid_annotation.fragment_newtork_dict.connect(self.on_fragment_network_loaded)

        self.data_export.create_menu_actions(self)
        self.data_export_menu.addAction(self.data_export.data_export_act)

        self.visualization.create_menu_actions(self)
        # TODO

    def init_widget(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget) 
        
        main_layout = QHBoxLayout(central_widget)
        
        self.data_widget = QTabWidget()
        self.left_layout = QVBoxLayout()
        
        raw_feature_widget = QWidget()
        raw_feature_layout = QVBoxLayout(raw_feature_widget)
        self.raw_feature_line_edit = QLineEdit()
        self.raw_feature_line_edit.setReadOnly(True)
        self.raw_feature_line_edit.mouseDoubleClickEvent = lambda event: self.on_raw_feature_double_click()
        raw_feature_layout.addWidget(self.raw_feature_line_edit)
        raw_feature_widget.setLayout(raw_feature_layout)
        
        annotated_feature_widget = QWidget()
        annotated_feature_layout = QVBoxLayout(annotated_feature_widget)
        self.annotated_feature_line_edit = QLineEdit()
        self.annotated_feature_line_edit.setReadOnly(True)
        self.annotated_feature_line_edit.mouseDoubleClickEvent = lambda event: self.on_annotated_feature_double_click()
        annotated_feature_layout.addWidget(self.annotated_feature_line_edit)
        annotated_feature_widget.setLayout(annotated_feature_layout)
        
        self.data_widget.addTab(raw_feature_widget, "Raw feature table")
        self.data_widget.addTab(annotated_feature_widget, "Annotated feature table")

        self.right_layout = QVBoxLayout()
    
        self.middle_widget = QTabWidget()
        self.middle_widget.setTabsClosable(True)
        self.middle_widget.tabCloseRequested.connect(self.close_tab)
        
        self.init_middle_widget(name="default table")
        
        self.right_widget = QWidget()
        
        self.init_right_widget()
        
        self.left_layout.addWidget(self.data_widget)
        self.right_layout.addWidget(self.middle_widget)
        self.right_layout.addWidget(self.right_widget)
        
        main_layout.addLayout(self.left_layout, 1)
        main_layout.addLayout(self.right_layout, 4)
        
        #self.centralWidget().setLayout(main_layout)
    
    def init_middle_widget(
        self,
        name: str,
    ):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        table = QTableWidget(10, 10)
        layout.addWidget(table)
        widget.setLayout(layout)

        self.middle_widget.addTab(widget, name)
        self.middle_widget.setCurrentIndex(self.middle_widget.count() - 1)

    def init_right_widget(self):
        """"""
        figure = Figure()
        canvas = FigureCanvas(figure)
        
        ax = figure.add_subplot(111)
        ax = none_spectrum_plot(
                title="",
                ax=ax
            )
        canvas.draw()
        self.right_widget = canvas

    def close_tab(
        self,
        index: int
    ):
        self.middle_widget.removeTab(index)
    
    def on_raw_ms1_peak_table_loaded(
        self,
        peak_table: pd.DataFrame
    ) -> None:
        self._ms1_peak_table = peak_table
        self.ms1_peak_table.emit(peak_table)
    
    def on_ms2_spectra_loaded(
        self,
        spectra: List[Spectrum]
    ) -> None:
        self._ms2_spectra = spectra
        self.ms2_spectra.emit(spectra)

    def on_reference_database_loaded(
        self,
        database: List[ReferenceLipid]
    ) -> None:
        self._reference_database = database
        self.reference_database.emit(database)
    
    def on_project_name_loaded(
        self,
        name: str
    ) -> None:
        self._project_name = name
        self.raw_feature_line_edit.setText(name)
    
    def on_annotated_peak_table_loaded(
        self,
        peak_table: pd.DataFrame
    ) -> None:
       self._peak_table_annotated = peak_table
       self.peak_table_annotated.emit(peak_table)
    
    def on_fragment_network_loaded(
        self,
        network_dict: Dict[str, FragmentNetowrk]
    ) -> None:
        self._fragment_newtork_dict = network_dict
        self.fragment_newtork_dict.emit(network_dict)
    
    def on_raw_feature_double_click(self) -> None:
        """"""    
        name = self._project_name + " raw feature table"
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        self.create_raw_feature_table()

        layout.addWidget(self.raw_table_widget)
        widget.setLayout(layout)

        self.middle_widget.addTab(widget, name)
        self.middle_widget.setCurrentIndex(self.middle_widget.count() - 1)
    
    def create_raw_feature_table(self) -> None:
        """"""
        if not hasattr(self, "_ms1_peak_table"):
            raise ValueError("No `ms1 peak table`")
        raw_data = self._ms1_peak_table.copy()
        num_row, num_col = raw_data.shape
        self.raw_table_widget = QTableWidget(num_row, num_col)
        
        items = [
            [QTableWidgetItem(str(val)) 
             for val in row] 
            for row in raw_data.values
        ]
        for row_idx, row_items in enumerate(items):
            for col_idx, item in enumerate(row_items):
                self.raw_table_widget.setItem(row_idx, col_idx, item)

        self.raw_table_widget.setHorizontalHeaderLabels(raw_data.columns.to_list())
        self.raw_table_widget.setSortingEnabled(True)
        self.raw_table_widget.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.raw_table_widget.cellDoubleClicked.connect(self.on_raw_peak_table_row_double_clicked)
    
    def on_raw_peak_table_row_double_clicked(
        self,
        row: int
    ) -> None:
        """"""
        data = self._ms1_peak_table.copy()
        feature_id_col = data.columns.get_loc("feature_id")
        feature_id = self.raw_table_widget.item(row, feature_id_col).text()
        primary_annotated_name_col = data.columns.get_loc("primary_annotated_name")
        primary_annotated_name = self.raw_table_widget.item(row, primary_annotated_name_col).text()
        title = f"MS/MS Spectrum of Feature {feature_id}"
        if primary_annotated_name:
            title = title + f"({primary_annotated_name})"
        
        spectrum_list = [
            spec
            for spec in self._ms2_spectra
            if (spec.metadata.get("id") == feature_id) or
            (spec.metadata.get("feature_id") == feature_id)
        ]
        
        figure = Figure()
        canvas = FigureCanvas(figure)
        ax = figure.add_subplot(111)
        if spectrum_list:
            spectrum = spectrum_list[0]
            ax = spectrum_plot(
                    spectrum=spectrum,
                    title=title,
                    grid=False,
                    ax=ax
                )
        else:
            ax = none_spectrum_plot(
                    title=title,
                    ax=ax
                )

        canvas.draw()

        if isinstance(self.right_widget, FigureCanvas):
            self.right_layout.removeWidget(self.right_widget)
            self.right_widget.deleteLater()

        self.right_layout.addWidget(canvas) 
        self.right_widget = canvas   
            
    def on_annotated_feature_double_click(self) -> None:
        """"""
        name = self._project_name + " annotated feature table"
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        self.create_annotated_feature_table()

        layout.addWidget(self.annotated_table_widget)
        widget.setLayout(layout)

        self.middle_widget.addTab(widget, name)
        self.middle_widget.setCurrentIndex(self.middle_widget.count() - 1)
    
    def create_annotated_feature_table(self) -> None:
        """
        """
        if not hasattr(self, "_peak_table_annotated"):
            raise ValueError("No `peak table annotated`")
    
        annotated_data = self._peak_table_annotated.copy()
        num_row, num_col = annotated_data.shape
        self.annotated_table_widget = QTableWidget(num_row, num_col)
        
        items = [
            [QTableWidgetItem(str(val)) 
             for val in row] 
            for row in annotated_data.values
        ]
        for row_idx, row_items in enumerate(items):
            for col_idx, item in enumerate(row_items):
                self.annotated_table_widget.setItem(row_idx, col_idx, item)

        self.annotated_table_widget.setHorizontalHeaderLabels(annotated_data.columns.to_list())
        self.annotated_table_widget.setSortingEnabled(True)
        self.annotated_table_widget.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.annotated_table_widget.cellDoubleClicked.connect(self.on_annotated_peak_table_row_double_clicked)

    def on_annotated_peak_table_row_double_clicked(
        self,
        row: int
    ):
        data = self._peak_table_annotated.copy()
        feature_id_col = data.columns.get_loc("feature_id")
        feature_id = self.raw_table_widget.item(row, feature_id_col).text()
        structure_annotated_name_col = data.columns.get_loc("structure_annotated_name")
        structure_annotated_name = self.raw_table_widget.item(row, structure_annotated_name_col).text()
        title = f"Fragment network of Feature {feature_id}"
        if structure_annotated_name:
            title = title + f"({structure_annotated_name})"

        fragment_network = self._fragment_newtork_dict.get(feature_id)
        
        figure = Figure()
        canvas = FigureCanvas(figure)
        ax = figure.add_subplot(111)
        ax = fragment_network_plot(
                fragment_network=fragment_network,
                structure_annotated_name=structure_annotated_name,
                ax=ax
            )
        canvas.draw()

        if isinstance(self.right_widget, FigureCanvas):
            self.right_layout.removeWidget(self.right_widget)
            self.right_widget.deleteLater()

        self.right_layout.addWidget(canvas) 
        self.right_widget = canvas   
        
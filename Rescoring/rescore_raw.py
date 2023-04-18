# ssh cadams@10.152.135.57
# conda activate oktoberfest-0_1
# python

import argparse
import logging
import os
import sys
import pandas as pd

from oktoberfest import __copyright__, __version__, logger, runner

import spectrum_fundamentals.constants as c
from spectrum_fundamentals.fragments import compute_peptide_mass
from spectrum_fundamentals.mod_string import internal_without_mods, maxquant_to_internal
from spectrum_io import Spectronaut
from spectrum_io.spectral_library import MSP

from oktoberfest.ce_calibration import CeCalibration, SpectralLibrary
from oktoberfest.data.spectra import Spectra
from oktoberfest.re_score import ReScore
from oktoberfest.utils.config import Config

from typing import List, Optional
from oktoberfest.utils.process_step import ProcessStep
from spectrum_io.raw import ThermoRaw
from oktoberfest.utils.multiprocessing_pool import JobPool
from spectrum_fundamentals.annotation.annotation import annotate_spectra
from oktoberfest.calculate_features import CalculateFeatures

search_dir = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/search-rescore-3"
config_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/search-rescore-3/test.json"

# runner.run_job(search_dir, config_path)

conf = Config()
conf.read(config_path)
job_type = conf.job_type
if conf.search_path:
    msms_path = conf.search_path
else:
    msms_path = os.path.join(search_dir, "msms.txt")

# run_rescoring(msms_path, search_dir, config_path)
re_score = ReScore(search_path=msms_path, raw_path=search_dir, out_path=search_dir, config_path=config_path)
re_score.get_raw_files()
re_score.split_msms()
# re_score.calculate_features()

# re_score.merge_input("rescore")
# re_score.merge_input("original")

# re_score.rescore_with_perc("rescore")
# re_score.rescore_with_perc("original")

mzml_path = re_score.get_mzml_folder_path()
# perc_path = re_score.get_percolator_folder_path()

for raw_file in re_score.raw_files:
    calc_feature_step = ProcessStep(re_score.out_path, "calculate_features." + raw_file)
    if calc_feature_step.is_done():
        continue

raw_file_path = os.path.join(re_score.raw_path, raw_file)
mzml_file_path = os.path.join(mzml_path, os.path.splitext(raw_file)[0] + ".mzML")
# percolator_input_path = re_score._get_split_perc_input_path(raw_file, "rescore")
split_msms_path = re_score._get_split_msms_path(raw_file)
# calculate_features_single(
#     raw_file_path,
#     split_msms_path,
#     percolator_input_path,
#     mzml_file_path,
#     config_path,
#     calc_feature_step,
# )
# if num_threads > 1:
#     processing_pool.check_pool(print_progress_every=1)


# def calculate_features_single
features = CalculateFeatures(search_path="", raw_path=raw_file_path, out_path=mzml_path, config_path=config_path)
df_search = pd.read_csv(split_msms_path, delimiter="\t")
# # features.predict_with_aligned_ce(df_search)

# def predict_with_aligned_ce
# # features.perform_alignment(df_search)

# def perform_alignment(self, df_search: pd.DataFrame):
hdf5_path = features.get_hdf5_path()
# if os.path.isfile(hdf5_path):
#     print(True)
# else:
#     print(False)
# # features.gen_lib(df_search)

# def gen_lib(self, df_search: Optional[pd.DataFrame] = None):
if df_search is None:
    raise AssertionError("You need to provide a dataframe.")
# df_raw = features._load_rawfile()

# # def _load_rawfile(self):

switch = features.config.raw_type
search_engine = features.config.search_type
logger.info(f"raw_type is {switch}")
# features._gen_mzml_from_thermo()
# def _gen_mzml_from_thermo(self):
raw = ThermoRaw()
# if not (features.out_path.endswith(".mzml")) and (not (features.out_path.endswith(".raw"))):
#     features.out_path = os.path.join(features.out_path, features.raw_path.split("/")[-1].split(".")[0] + ".mzml")
# features.raw_path = raw.convert_raw_mzml(input_path=features.raw_path, output_path=features.out_path)
# features.out_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/search-rescore-3/mzML/220329_NHG_benign_UDN31_PBMC_W6-32_17%_orbitrap_DDA_Rep2.mzml"
print(features.raw_path)
# features.raw_path = features.raw_path.as_posix().replace(".raw", ".mzml")
df_raw = ThermoRaw.read_mzml(source=features.out_path, package=features.mzml_reader_package, search_type=search_engine)

df_raw_drop = df_raw.drop(columns = ['MASS_ANALYZER', 'FRAGMENTATION'])
df_search[['MASS_ANALYZER', 'FRAGMENTATION']]
df_raw[['MASS_ANALYZER', 'FRAGMENTATION']]

df_join = df_search.merge(df_raw_drop, on=["RAW_FILE", "SCAN_NUMBER"])
logger.info(f"There are {len(df_join)} matched identifications")
df_join.columns
logger.info("Annotating raw spectra")
df_annotated_spectra = annotate_spectra(df_join)

df_join.drop(columns=["INTENSITIES", "MZ"], inplace=True)
# return df_annotated_spectra["INTENSITIES"]
logger.info("Preparing library")
self.library.add_columns(df_join)
self.library.add_matrix(df_annotated_spectra["INTENSITIES"], FragmentType.RAW)
self.library.add_matrix(df_annotated_spectra["MZ"], FragmentType.MZ)
self.library.add_column(df_annotated_spectra["CALCULATED_MASS"], "CALCULATED_MASS")


# return df_search
logger.info("Merging rawfile and search result")
df_join = df_search.merge(df_raw, on=["RAW_FILE", "SCAN_NUMBER"])
logger.info(f"There are {len(df_join)} matched identifications")
logger.info("Annotating raw spectra")
df_annotated_spectra = annotate_spectra(df_join)
df_join.drop(columns=["INTENSITIES", "MZ"], inplace=True)
# return df_annotated_spectra["INTENSITIES"]
logger.info("Preparing library")
self.library.add_columns(df_join)
self.library.add_matrix(df_annotated_spectra["INTENSITIES"], FragmentType.RAW)
self.library.add_matrix(df_annotated_spectra["MZ"], FragmentType.MZ)
self.library.add_column(df_annotated_spectra["CALCULATED_MASS"], "CALCULATED_MASS")





features.library.spectra_data["COLLISION_ENERGY"] = features.best_ce
features.grpc_predict(features.library)
features.library.write_pred_as_hdf5(features.get_pred_path())


hdf5_path = features.get_hdf5_path()
if os.path.isfile(hdf5_path):
    print(True)
    features.library.read_from_hdf5(hdf5_path)
else:
    features.gen_lib(df_search)
    features.write_metadata_annotation()

df_raw = features._load_rawfile()

switch = features.config.raw_type
search_engine = features.config.search_type
if switch == "thermo":
    features._gen_mzml_from_thermo()








# return df_search
logger.info("Merging rawfile and search result")
df_join = df_search.merge(df_raw, on=["RAW_FILE", "SCAN_NUMBER"])
logger.info(f"There are {len(df_join)} matched identifications")

logger.info("Annotating raw spectra")
df_annotated_spectra = annotate_spectra(df_join)
df_join.drop(columns=["INTENSITIES", "MZ"], inplace=True)
# return df_annotated_spectra["INTENSITIES"]
logger.info("Preparing library")
self.library.add_columns(df_join)
self.library.add_matrix(df_annotated_spectra["INTENSITIES"], FragmentType.RAW)
self.library.add_matrix(df_annotated_spectra["MZ"], FragmentType.MZ)
self.library.add_column(df_annotated_spectra["CALCULATED_MASS"], "CALCULATED_MASS")


logger = logging.getLogger(__name__)

def generate_spectral_lib(search_dir: str, config_path: str):
    """
    Create a SpectralLibrary object and generate the spectral library.

    :param search_dir: path to directory containing the msms.txt and raw files
    :param config_path: path to config file
    :raises ValueError: spectral library output format is not supported as spectral library type
    """
    spec_library = SpectralLibrary(path=search_dir, out_path=search_dir, config_path=config_path)
    spec_library.gen_lib()
    spec_library.library.spectra_data["MODIFIED_SEQUENCE"] = spec_library.library.spectra_data[
        "MODIFIED_SEQUENCE"
    ].apply(lambda x: "_" + x + "_")
    models_dict = spec_library.config.models
    spec_library.library.spectra_data["MODIFIED_SEQUENCE"] = maxquant_to_internal(
        spec_library.library.spectra_data["MODIFIED_SEQUENCE"], fixed_mods={}
    )
    spec_library.library.spectra_data["SEQUENCE"] = internal_without_mods(
        spec_library.library.spectra_data["MODIFIED_SEQUENCE"]
    )
    spec_library.library.spectra_data["PEPTIDE_LENGTH"] = spec_library.library.spectra_data["SEQUENCE"].apply(
        lambda x: len(x)
    )

    logger.info(f"No of sequences before Filtering is {len(spec_library.library.spectra_data['PEPTIDE_LENGTH'])}")
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        (spec_library.library.spectra_data["PEPTIDE_LENGTH"] <= 30)
    ]
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        (~spec_library.library.spectra_data["MODIFIED_SEQUENCE"].str.contains(r"\(ac\)"))
    ]
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        (~spec_library.library.spectra_data["MODIFIED_SEQUENCE"].str.contains(r"\(Acetyl \(Protein N-term\)\)"))
    ]
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        (~spec_library.library.spectra_data["SEQUENCE"].str.contains("U"))
    ]
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        spec_library.library.spectra_data["PRECURSOR_CHARGE"] <= 6
    ]
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        spec_library.library.spectra_data["PEPTIDE_LENGTH"] >= 7
    ]
    logger.info(f"No of sequences after Filtering is {len(spec_library.library.spectra_data['PEPTIDE_LENGTH'])}")

    tmt_model = False
    for _, value in models_dict.items():
        if value:
            if "TMT" in value:
                tmt_model = True
    if tmt_model and spec_library.config.tag != "":
        unimod_tag = c.TMT_MODS[spec_library.config.tag]
        spec_library.library.spectra_data["MODIFIED_SEQUENCE"] = maxquant_to_internal(
            spec_library.library.spectra_data["MODIFIED_SEQUENCE"],
            fixed_mods={"C": "C[UNIMOD:4]", "^_": f"_{unimod_tag}", "K": f"K{unimod_tag}"},
        )
    else:
        spec_library.library.spectra_data["MODIFIED_SEQUENCE"] = maxquant_to_internal(
            spec_library.library.spectra_data["MODIFIED_SEQUENCE"]
        )
    spec_library.library.spectra_data["MASS"] = spec_library.library.spectra_data["MODIFIED_SEQUENCE"].apply(
        lambda x: compute_peptide_mass(x)
    )
    no_of_spectra = len(spec_library.library.spectra_data)
    no_of_sections = no_of_spectra // 7000
    for i in range(0, no_of_sections + 1):
        spectra_div = Spectra()
        if i < no_of_sections:
            spectra_div.spectra_data = spec_library.library.spectra_data.iloc[i * 7000 : (i + 1) * 7000]
            logger.info(f"Indices {i * 7000}, {(i + 1) * 7000}")
        elif (i * 7000) < no_of_spectra:
            spectra_div.spectra_data = spec_library.library.spectra_data.iloc[i * 7000 :]
            logger.info(f"Last Batch from index {i * 7000}")
            logger.info(f"Batch of size {len(spectra_div.spectra_data.index)}")
        else:
            break

        grpc_output_sec = spec_library.grpc_predict(spectra_div)
        if spec_library.config.output_format == "msp":
            out_lib_msp = MSP(
                spectra_div.spectra_data, grpc_output_sec, os.path.join(spec_library.results_path, "myPrositLib.msp")
            )
            out_lib_msp.prepare_spectrum()
            out_lib_msp.write()
        elif spec_library.config.output_format == "spectronaut":
            out_lib_spectronaut = Spectronaut(
                spectra_div.spectra_data, grpc_output_sec, os.path.join(spec_library.results_path, "myPrositLib.csv")
            )
            out_lib_spectronaut.prepare_spectrum()
            out_lib_spectronaut.write()
        else:
            raise ValueError(f"{spec_library.config.output_format} is not supported as spectral library type")


def run_ce_calibration(msms_path: str, search_dir: str, config_path: str):
    """
    Create a CeCalibration object and run the CE calibration.

    :param msms_path: path to msms folder
    :param search_dir: path to directory containing the msms.txt and raw files
    :param config_path: path to config file
    :raises ValueError: raw_type is not supported as rawfile-type
    """
    ce_calib = CeCalibration(search_path=msms_path, raw_path=search_dir, out_path=search_dir, config_path=config_path)
    df_search = ce_calib._load_search()
    raw_type = ce_calib.config.raw_type
    if raw_type == "thermo":
        extension = ".raw"
    elif raw_type == "mzml":
        extension = ".mzml"
    else:
        raise ValueError(f"{raw_type} is not supported as rawfile-type")

    ce_calib.raw_path = os.path.join(
        ce_calib.raw_path,
        [os.path.basename(f) for f in os.listdir(ce_calib.raw_path) if f.lower().endswith(extension)][0],
    )
    ce_calib.perform_alignment(df_search)
    with open(os.path.join(ce_calib.results_path, "ce.txt"), "w") as f:
        f.write(str(ce_calib.best_ce))




class CalculateFeatures(CeCalibration):
    """
    Main to init a re-score obj and go through the steps.
    1- predict_with_aligned_ce
    2- gen_perc_metrics
    """
    def predict_with_aligned_ce(self, df_search: pd.DataFrame):
        """
        Get best collision energy with ce_calibration then use it for prediction.
        :param df_search: a msms matrix as a pd.DataFrame
        """
        self.perform_alignment(df_search)
        self.library.spectra_data["COLLISION_ENERGY"] = self.best_ce
        self.grpc_predict(self.library)
        self.library.write_pred_as_hdf5(self.get_pred_path())
    def gen_perc_metrics(self, search_type: str, file_path: Optional[str]):
        """
        Get all percolator metrics and add them to library.
        :param search_type: model (rescore or original) as a string
        :param file_path: path to percolator input file as a string
        """
        perc_features = Percolator(
            metadata=self.library.get_meta_data(),
            pred_intensities=self.library.get_matrix(FragmentType.PRED),
            true_intensities=self.library.get_matrix(FragmentType.RAW),
            input_type=search_type,
            all_features_flag=self.config.all_features,
        )
        perc_features.calc()
        if file_path:
            perc_features.write_to_file(file_path)

dir(CeCalibration)

class CeCalibration(SpectralLibrary):
    """
    Main to init a CeCalibrarion obj and go through the steps.
    1- gen_lib
    2- allign_ce
    3- get_best_ce
    4- write output
    """
    raw_path: str
    out_path: str
    best_ce: float
    def __init__(
        self,
        search_path: str,
        raw_path: str,
        out_path: str,
        config_path: Optional[str],
        mzml_reader_package: str = "pyteomics",
    ):
        """
        Initialize a CeCalibration object.
        :param search_path: path to search directory
        :param raw_path: path to directory containing the msms.txt and raw files
        :param out_path: path to output folder
        :param config_path: path to configuration file
        :param mzml_reader_package: mzml reader (pymzml or pyteomics)
        """
        super().__init__(search_path, out_path, config_path=config_path)
        self.search_path = search_path
        self.raw_path = raw_path
        self.out_path = out_path
        self.mzml_reader_package = mzml_reader_package
        self.best_ce = 0
    def _gen_internal_search_result_from_msms(self):
        """Generate internal search result from msms.txt."""
        logger.info(f"Converting msms.txt at location {self.search_path} to internal search result.")
        models_dict = self.config.models
        tmt_model = False
        for _, value in models_dict.items():
            if value and "TMT" in value:
                tmt_model = True
        if tmt_model:
            tmt_labeled = self.config.tag
        else:
            tmt_labeled = ""
        search_type = self.config.search_type
        if search_type == "Maxquant":
            mxq = MaxQuant(self.search_path)
            self.search_path = mxq.generate_internal(tmt_labeled=tmt_labeled)
        elif search_type == "Msfragger":
            msf = MSFragger(self.search_path)
            self.search_path = msf.generate_internal(tmt_labeled=tmt_labeled)
        elif search_type == "Mascot":
            mascot = Mascot(self.search_path)
            self.search_path = mascot.generate_internal(tmt_labeled=tmt_labeled)
    def _gen_mzml_from_thermo(self):
        """Generate mzml from thermo raw file."""
        logger.info("Converting thermo rawfile to mzml.")
        raw = ThermoRaw()
        if not (self.out_path.endswith(".mzML")) and (not (self.out_path.endswith(".raw"))):
            self.out_path = os.path.join(self.out_path, self.raw_path.split("/")[-1].split(".")[0] + ".mzml")
        self.raw_path = raw.convert_raw_mzml(input_path=self.raw_path, output_path=self.out_path)
    def _load_search(self):
        """Load search type."""
        switch = self.config.search_type
        logger.info(f"search_type is {switch}")
        if switch == "Maxquant" or switch == "Msfragger" or switch == "Mascot":
            self._gen_internal_search_result_from_msms()
        elif switch == "Internal":
            pass
        else:
            raise ValueError(f"{switch} is not supported as search-type")
        if switch == "Maxquant":
            return MaxQuant.read_internal(MaxQuant(self.search_path), self.search_path)
        elif switch == "Msfragger":
            return MSFragger.read_internal(MSFragger(self.search_path), path=self.search_path)
        else:
            return Mascot.read_internal(Mascot(self.search_path), path=self.search_path)
    def _load_rawfile(self):
        """Load raw file."""
        switch = self.config.raw_type
        search_engine = self.config.search_type
        logger.info(f"raw_type is {switch}")
        if switch == "thermo":
            self._gen_mzml_from_thermo()
        elif switch == "mzml":
            pass
        else:
            raise ValueError(f"{switch} is not supported as rawfile-type")
        print(self.raw_path)
        self.raw_path = self.raw_path.as_posix().replace(".raw", ".mzml")
        return ThermoRaw.read_mzml(source=self.out_path, package=self.mzml_reader_package, search_type=search_engine)
    def gen_lib(self, df_search: Optional[pd.DataFrame] = None):
        """
        Read input search and raw and add it to library.
        Method inherits from superclass, therefore the Optional is required to ensure same method signature.
        It needs to be refactored in the future.
        :param df_search: search result as pd.DataFrame
        :raises AssertionError: raises if df_search is not given
        """
        if df_search is None:
            raise AssertionError("You need to provide a dataframe.")
        df_raw = self._load_rawfile()
        # return df_search
        logger.info("Merging rawfile and search result")
        df_join = df_search.merge(df_raw, on=["RAW_FILE", "SCAN_NUMBER"])
        logger.info(f"There are {len(df_join)} matched identifications")
        logger.info("Annotating raw spectra")
        df_annotated_spectra = annotate_spectra(df_join)
        df_join.drop(columns=["INTENSITIES", "MZ"], inplace=True)
        # return df_annotated_spectra["INTENSITIES"]
        logger.info("Preparing library")
        self.library.add_columns(df_join)
        self.library.add_matrix(df_annotated_spectra["INTENSITIES"], FragmentType.RAW)
        self.library.add_matrix(df_annotated_spectra["MZ"], FragmentType.MZ)
        self.library.add_column(df_annotated_spectra["CALCULATED_MASS"], "CALCULATED_MASS")
    def get_hdf5_path(self) -> str:
        """Get path to hdf5 file."""
        return self.out_path + ".hdf5"
    def get_pred_path(self) -> str:
        """Get path to prediction hdf5 file."""
        return self.out_path + "_pred.hdf5"
    def write_metadata_annotation(self):
        """Write metadata annotation as hdf5 file."""
        self.library.write_as_hdf5(self.get_hdf5_path())
    def _prepare_alignment_df(self):
        self.alignment_library = Spectra()
        self.alignment_library.spectra_data = self.library.spectra_data.copy()
        # Remove decoy and HCD fragmented spectra
        self.alignment_library.spectra_data = self.alignment_library.spectra_data[
            (self.alignment_library.spectra_data["FRAGMENTATION"] == "HCD")
            & (~self.alignment_library.spectra_data["REVERSE"])
        ]
        # Select the 1000 highest scoring or all if there are less than 1000
        self.alignment_library.spectra_data = self.alignment_library.spectra_data.sort_values(
            by="SCORE", ascending=False
        ).iloc[:1000]
        # Repeat dataframe for each CE
        ce_range = range(18, 50)
        nrow = len(self.alignment_library.spectra_data)
        self.alignment_library.spectra_data = pd.concat([self.alignment_library.spectra_data for _ in ce_range], axis=0)
        self.alignment_library.spectra_data["COLLISION_ENERGY"] = np.repeat(ce_range, nrow)
        self.alignment_library.spectra_data.reset_index(inplace=True)
    def _predict_alignment(self):
        self.grpc_predict(self.alignment_library, alignment=True)
    def _alignment(self):
        """
        Edit library to try different ranges of ce 15-50. then predict with the new library.
        Check https://gitlab.lrz.de/proteomics/prosit_tools/oktoberfest/-/blob/develop/oktoberfest/ce_calibration/grpc_alignment.py
        """
        pred_intensity = self.alignment_library.get_matrix(FragmentType.PRED)
        raw_intensity = self.alignment_library.get_matrix(FragmentType.RAW)
        # return pred_intensity.toarray(), raw_intensity.toarray()
        sm = SimilarityMetrics(pred_intensity, raw_intensity)
        self.alignment_library.spectra_data["SPECTRAL_ANGLE"] = sm.spectral_angle(raw_intensity, pred_intensity, 0)
        self.ce_alignment = self.alignment_library.spectra_data.groupby(by=["COLLISION_ENERGY"])[
            "SPECTRAL_ANGLE"
        ].mean()
        if "/" in self.raw_path:
            split_char = "/"
        else:
            split_char = "\\"
        plot_mean_sa_ce(
            sa_ce_df=self.ce_alignment,
            directory=os.path.join((split_char).join(self.raw_path.split(split_char)[:-1]), "results"),
            raw_file_name=self.raw_path.split(split_char)[-1],
        )
    def _get_best_ce(self):
        """Get aligned ce for this lib."""
        self.best_ce = self.ce_alignment.idxmax()
        logger.info(f"Best collision energy: {self.best_ce}")
    def perform_alignment(self, df_search: pd.DataFrame):
        """
        Perform alignment and get the best CE.
        :param df_search: search result as pd.DataFrame
        """
        hdf5_path = self.get_hdf5_path()
        logger.info(f"Path to hdf5 file with annotations for {self.out_path}: {hdf5_path}")
        if os.path.isfile(hdf5_path):
            self.library.read_from_hdf5(hdf5_path)
        else:
            self.gen_lib(df_search)
            self.write_metadata_annotation()
        # Check if all data is HCD no need to align and return the best ce as 35
        hcd_df = self.library.spectra_data[(self.library.spectra_data["FRAGMENTATION"] == "HCD")]
        if len(hcd_df.index) == 0:
            self.best_ce = 35.0
            return
        self._prepare_alignment_df()
        self._predict_alignment()
        self._alignment()
        self._get_best_ce()


class SpectralLibrary:
    """
    Main to init a SpectralLibrary obj and go through the steps.
    1- gen_lib
    2- grpc_predict
    3- write output
    """
    path: str
    library: Spectra
    config: Config
    config_path: Optional[str]
    num_threads: int
    grpc_output: dict
    def __init__(self, path: str, out_path: str, config_path: Optional[str]):
        """
        Initialize a SpectralLibrary object.
        :param path: path to directory containing the msms.txt and raw files
        :param out_path: path to output folder
        :param config_path: path to configuration file
        """
        self.path = path
        self.library = Spectra()
        self.config_path = config_path
        self.config = Config()
        if config_path:
            self.config.read(config_path)
        else:
            self.config.read(CONFIG_PATH)
        self.results_path = os.path.join(out_path, "results")
        if os.path.isdir(out_path):
            if not os.path.isdir(self.results_path):
                try:
                    os.makedirs(self.results_path)
                except Exception:
                    print("In Feature Calculation")
        else:
            print("In Feature Calculation")
    def gen_lib(self, df_search: Optional[pd.DataFrame] = None):
        """
        Read input csv file and add it to library.
        :param df_search: unused, necessary to ensure same method signature for inheriting function
        """
        if self.config.fasta:
            self.read_fasta()
            library_df = csv.read_file(os.path.join(self.path, "prosit_input.csv"))
        else:
            for file in os.listdir(self.path):
                if file.endswith(".csv"):
                    library_df = csv.read_file(os.path.join(self.path, file))
        library_df.columns = library_df.columns.str.upper()
        self.library.add_columns(library_df)
    def grpc_predict(self, library: Spectra, alignment: bool = False):
        """
        Use grpc to predict library and add predictions to library.
        :param library: Spectra object with the library
        :param alignment: True if alignment present
        :return: grpc predictions if we are trying to generate spectral library
        """
        path = Path(__file__).parent / "certificates/"
        logger.info(path)
        predictor = PROSITpredictor(
            server=self.config.prosit_server,
            path_to_ca_certificate=os.path.join(path, "Proteomicsdb-Prosit-v2.crt"),
            path_to_certificate=os.path.join(path, "oktoberfest-production.crt"),
            path_to_key_certificate=os.path.join(path, "oktoberfest-production.key"),
        )
        models_dict = self.config.models
        models = []
        tmt_model = False
        for _, value in models_dict.items():
            if not value:
                continue
            tmt_model = True if "TMT" in value else tmt_model
            models.append(value)
            if alignment:
                break
        if tmt_model:
            library.spectra_data["FRAGMENTATION_GRPC"] = library.spectra_data["FRAGMENTATION"].apply(
                lambda x: 2 if x == "HCD" else 1
            )
        library.spectra_data["GRPC_SEQUENCE"] = library.spectra_data["MODIFIED_SEQUENCE"]
        try:
            predictions = predictor.predict(
                sequences=library.spectra_data["GRPC_SEQUENCE"].values.tolist(),
                charges=library.spectra_data["PRECURSOR_CHARGE"].values.tolist(),
                collision_energies=library.spectra_data["COLLISION_ENERGY"].values / 100.0,
                fragmentation=library.spectra_data["FRAGMENTATION_GRPC"].values if tmt_model else None,
                models=models,
                disable_progress_bar=True,
            )
        except BaseException:
            logger.exception("An exception was thrown!", exc_info=True)
            print(library.spectra_data["GRPC_SEQUENCE"])
        # Return only in spectral library generation otherwise add to library
        if self.config.job_type == "SpectralLibraryGeneration":
            return predictions
        intensities_pred = pd.DataFrame()
        intensities_pred["intensity"] = predictions[models[0]]["intensity"].tolist()
        library.add_matrix(intensities_pred["intensity"], FragmentType.PRED)
        if alignment:
            return
        irt_pred = predictions[models[1]]
        library.add_column(irt_pred, "PREDICTED_IRT")
        if len(models) > 2:
            proteotypicity_pred = predictions[models[2]]
            library.add_column(proteotypicity_pred, "PROTEOTYPICITY")
    def read_fasta(self):
        """Read fasta file."""
        cmd = [
            "--fasta",
            f"{self.config.fasta}",
            "--prosit_input",
            f"{os.path.join(self.path, 'prosit_input.csv')}",
            "--fragmentation",
            f"{self.config.fragmentation}",
            "--digestion",
            f"{self.config.digestion}",
            "--cleavages",
            f"{self.config.cleavages}",
            "--db",
            f"{self.config.db}",
            "--enzyme",
            f"{self.config.enzyme}",
            "--special-aas",
            f"{self.config.special_aas}",
            "--min-length",
            f"{self.config.min_length}",
            "--max-length",
            f"{self.config.max_length}",
        ]
        digest.main(cmd)



class ReScore(CalculateFeatures):
    """
    Main to init a re-score obj and go through the steps.
    1- get_raw_files
    2- split_msms
    3- calculate_features
    4- merge_input
    5- rescore_with_perc
    """
    raw_files: List[str]
    split_msms_step: ProcessStep
    merge_input_step_prosit: ProcessStep
    merge_input_step_andromeda: ProcessStep
    percolator_step_prosit: ProcessStep
    percolator_step_andromeda: ProcessStep
    def __init__(
        self,
        search_path: str,
        raw_path: str,
        out_path: str,
        config_path: Optional[str],
        mzml_reader_package: str = "pymzml",
    ):
        """
        Initialize a ReScore object and go through the steps.
        1- get_raw_files
        2- split_msms
        3- calculate_features
        4- merge_input
        5- rescore_with_perc
        :param search_path: path to search directory
        :param raw_path: path to raw file as a string
        :param out_path: path to output folder
        :param config_path: path to config file
        :param mzml_reader_package: mzml reader (pymzml or pyteomics)
        """
        super().__init__(
            search_path, raw_path, out_path, config_path=config_path, mzml_reader_package=mzml_reader_package
        )
        self.split_msms_step = ProcessStep(out_path, "split_msms")
        self.merge_input_step_prosit = ProcessStep(out_path, "merge_input_prosit")
        self.merge_input_step_andromeda = ProcessStep(out_path, "merge_input_andromeda")
        self.percolator_step_prosit = ProcessStep(out_path, "percolator_prosit")
        self.percolator_step_andromeda = ProcessStep(out_path, "percolator_andromeda")
    def get_raw_files(self):
        """
        Obtains raw files by scanning through the raw_path directory.
        If raw_path is a file, only process this one.
        :raises ValueError: raw_type is not supported as rawfile-type
        """
        self.raw_files = []
        if os.path.isfile(self.raw_path):
            self.raw_files = [self.raw_path]
            self.raw_path = os.path.dirname(self.raw_path)
        elif os.path.isdir(self.raw_path):
            raw_type = self.config.raw_type
            if raw_type == "thermo":
                extension = ".raw"
            elif raw_type == "mzml":
                extension = ".mzml"
            else:
                raise ValueError(f"{raw_type} is not supported as rawfile-type")
            self.raw_files = [os.path.basename(f) for f in os.listdir(self.raw_path) if f.lower().endswith(extension)]
            logger.info(f"Found {len(self.raw_files)} raw files in the search directory")
    def split_msms(self):
        """Splits msms.txt file per raw file such that we can process each raw file in parallel \
        without reading the entire msms.txt."""
        if self.split_msms_step.is_done():
            return
        msms_path = self.get_msms_folder_path()
        if not os.path.isdir(msms_path):
            os.makedirs(msms_path)
        df_search = self._load_search()
        logger.info(f"Read {len(df_search.index)} PSMs from {self.search_path}")
        for raw_file, df_search_split in df_search.groupby("RAW_FILE"):
            raw_file_path = os.path.join(self.raw_path, raw_file)
            if not (os.path.isfile(raw_file_path + ".raw") or os.path.isfile(raw_file_path + ".RAW")):
                logger.info(f"Did not find {raw_file} in search directory, skipping this file")
                continue
            split_msms = self._get_split_msms_path(raw_file)
            logger.info(f"Creating split msms.txt file {split_msms}")
            df_search_split = df_search_split[(df_search_split["PEPTIDE_LENGTH"] <= 30)]
            df_search_split = df_search_split[(~df_search_split["MODIFIED_SEQUENCE"].str.contains(r"\(ac\)"))]
            df_search_split = df_search_split[
                (~df_search_split["MODIFIED_SEQUENCE"].str.contains(r"\(Acetyl \(Protein N-term\)\)"))
            ]
            df_search_split = df_search_split[(~df_search_split["SEQUENCE"].str.contains("U"))]
            df_search_split = df_search_split[df_search_split["PRECURSOR_CHARGE"] <= 6]
            df_search_split = df_search_split[df_search_split["PEPTIDE_LENGTH"] >= 7]
            df_search_split.to_csv(split_msms, sep="\t", index=False)
        self.split_msms_step.mark_done()
    def calculate_features(self):
        """Calculates percolator input features per raw file using multiprocessing."""
        num_threads = self.config.num_threads
        self.config
        if num_threads > 1:
            processing_pool = JobPool(processes=num_threads)
        mzml_path = self.get_mzml_folder_path()
        if not os.path.isdir(mzml_path):
            os.makedirs(mzml_path)
        perc_path = self.get_percolator_folder_path()
        if not os.path.isdir(perc_path):
            os.makedirs(perc_path)
        for raw_file in self.raw_files:
            calc_feature_step = ProcessStep(self.out_path, "calculate_features." + raw_file)
            if calc_feature_step.is_done():
                continue
            raw_file_path = os.path.join(self.raw_path, raw_file)
            mzml_file_path = os.path.join(mzml_path, os.path.splitext(raw_file)[0] + ".mzML")
            percolator_input_path = self._get_split_perc_input_path(raw_file, "rescore")
            split_msms_path = self._get_split_msms_path(raw_file)
            if num_threads > 1:
                processing_pool.apply_async(
                    calculate_features_single,
                    (
                        raw_file_path,
                        split_msms_path,
                        percolator_input_path,
                        mzml_file_path,
                        self.config_path,
                        calc_feature_step,
                    ),
                )
            else:
                calculate_features_single(
                    raw_file_path,
                    split_msms_path,
                    percolator_input_path,
                    mzml_file_path,
                    self.config_path,
                    calc_feature_step,
                )
        if num_threads > 1:
            processing_pool.check_pool(print_progress_every=1)
    def merge_input(self, search_type: str = "rescore"):
        """
        Merge percolator input files into one large file for combined percolation.
        Fastest solution according to:
        https://stackoverflow.com/questions/44211461/what-is-the-fastest-way-to-combine-100-csv-files-with-headers-into-one
        :param search_type: choose either rescore or original to merge percolator files for this.
        """
        if search_type == "rescore":
            if self.merge_input_step_prosit.is_done():
                return
        else:
            if self.merge_input_step_andromeda.is_done():
                return
        merged_perc_input_file_prosit = self._get_merged_perc_input_path(search_type)
        logger.info("Merging percolator input files for " + search_type)
        with open(merged_perc_input_file_prosit, "wb") as fout:
            first = True
            for raw_file in self.raw_files:
                percolator_input_path = self._get_split_perc_input_path(raw_file, search_type)
                with open(percolator_input_path, "rb") as f:
                    if not first:
                        next(f)  # skip the header
                    else:
                        first = False
                    fout.write(f.read())
                os.remove(percolator_input_path)
        df_prosit = pd.read_csv(merged_perc_input_file_prosit, sep="\t")
        df_prosit = df_prosit.fillna(0)
        df_prosit.to_csv(merged_perc_input_file_prosit, sep="\t", index=False)
        if search_type == "rescore":
            self.merge_input_step_prosit.mark_done()
        else:
            self.merge_input_step_andromeda.mark_done()
    def rescore_with_perc(self, search_type: str = "rescore", test_fdr: float = 0.01, train_fdr: float = 0.01):
        """Use percolator to re-score library."""
        if search_type == "rescore":
            if self.percolator_step_prosit.is_done():
                return
        else:
            if self.percolator_step_andromeda.is_done():
                return
        perc_path = self.get_percolator_folder_path()
        weights_file = os.path.join(perc_path, f"{search_type}_weights.csv")
        target_psms = os.path.join(perc_path, f"{search_type}_target.psms")
        decoy_psms = os.path.join(perc_path, f"{search_type}_decoy.psms")
        target_peptides = os.path.join(perc_path, f"{search_type}_target.peptides")
        decoy_peptides = os.path.join(perc_path, f"{search_type}_decoy.peptides")
        log_file = os.path.join(perc_path, f"{search_type}.log")
        cmd = f"percolator --weights {weights_file} \
                          --num-threads {self.config.num_threads} \
                          --subset-max-train 500000 \
                          --post-processing-tdc \
                          --search-input concatenated \
                          --testFDR {test_fdr} \
                          --trainFDR {train_fdr} \
                          --results-psms {target_psms} \
                          --decoy-results-psms {decoy_psms} \
                          --results-peptides {target_peptides} \
                          --decoy-results-peptides {decoy_peptides} \
                          {self._get_merged_perc_input_path(search_type)} 2> {log_file}"
        logger.info(f"Starting percolator with command {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        if search_type == "rescore":
            self.percolator_step_prosit.mark_done()
        else:
            plot_all(perc_path)
            self.percolator_step_andromeda.mark_done()
    def get_msms_folder_path(self):
        """Get folder path to msms."""
        return os.path.join(self.out_path, "msms")
    def _get_split_msms_path(self, raw_file: str) -> str:
        """
        Get path to split msms.
        :param raw_file: path to raw file as a string
        :return: path to split msms file
        """
        return os.path.join(self.get_msms_folder_path(), os.path.splitext(raw_file)[0] + ".rescore")
    def get_mzml_folder_path(self) -> str:
        """Get folder path to mzml."""
        return os.path.join(self.out_path, "mzML")
    def get_percolator_folder_path(self) -> str:
        """Get folder path to percolator."""
        return os.path.join(self.results_path, "percolator")
    def _get_split_perc_input_path(self, raw_file: str, search_type: str):
        """
        Specify search_type to differentiate between percolator and andromeda output.
        :param raw_file: path to raw file as a string
        :param search_type: model (rescore or original) as a string
        :return: path to split percolator input file
        """
        return os.path.join(
            self.get_percolator_folder_path(), os.path.splitext(raw_file)[0] + "_" + search_type + ".tab"
        )
    def _get_merged_perc_input_path(self, search_type: str):
        """
        Get merged percolator input path.
        :param search_type: model (rescore or original) as a string
        :return: path to merged percolator input folder
        """
        return os.path.join(self.get_percolator_folder_path(), search_type + ".tab")


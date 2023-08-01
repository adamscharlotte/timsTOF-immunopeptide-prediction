# ssh cadams@10.152.135.57
# conda activate oktoberfest-0_1
# python

import logging
import os

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

from oktoberfest.ce_calibration import SpectralLibrary
from oktoberfest.utils.process_step import ProcessStep
from oktoberfest.calculate_features import CalculateFeatures
import pandas as pd
from prosit_grpc.predictPROSIT import PROSITpredictor
from pathlib import Path

from tensorflow_serving.apis import model_service_pb2_grpc
from tensorflow_serving.apis import get_model_status_pb2

search_dir = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/search-test-2"
config_path = "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/search-test-2/test.json"

runner.run_job(search_dir, config_path)


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

# --------- calculate_features --------
num_threads = re_score.config.num_threads
re_score.config
mzml_path = re_score.get_mzml_folder_path()
perc_path = re_score.get_percolator_folder_path()

for raw_file in re_score.raw_files:
    calc_feature_step = ProcessStep(re_score.out_path, "calculate_features." + raw_file)
    if calc_feature_step.is_done():
        continue

raw_file_path = os.path.join(re_score.raw_path, raw_file)
mzml_file_path = os.path.join(mzml_path, os.path.splitext(raw_file)[0] + ".mzML")
percolator_input_path = re_score._get_split_perc_input_path(raw_file, "rescore")
split_msms_path = re_score._get_split_msms_path(raw_file)
# calculate_features_single(
#     raw_file_path,
#     split_msms_path,
#     percolator_input_path,
#     mzml_file_path,
#     config_path,
#     calc_feature_step,
# )

# --------- calculate_features_single ---------
features = CalculateFeatures(search_path="", raw_path=raw_file_path, out_path=mzml_path, config_path=config_path)
df_search = pd.read_csv(split_msms_path, delimiter="\t")
features.predict_with_aligned_ce(df_search)

# # --------- predict_with_aligned_ce ---------
# features.perform_alignment(df_search)
# features.library.spectra_data["COLLISION_ENERGY"] = features.best_ce
# features.grpc_predict(features.library)
# features.library.write_pred_as_hdf5(features.get_pred_path())
features.gen_perc_metrics("rescore", percolator_input_path)

# --------- gen_perc_metrics ---------
from spectrum_fundamentals.metrics.percolator import Percolator
from oktoberfest.data.spectra import FragmentType

perc_features_test_2 = Percolator(
        metadata=features.library.get_meta_data(),
        pred_intensities=features.library.get_matrix(FragmentType.PRED),
        true_intensities=features.library.get_matrix(FragmentType.RAW),
        mz=features.library.get_matrix(FragmentType.MZ),
        input_type="rescore",
        all_features_flag=features.config.all_features,
    )
perc_features_test_2.calc()
dir(perc_features)

# --------- calc ---------
import numpy as np
from spectrum_fundamentals.metrics import fragments_ratio as fr
from spectrum_fundamentals.metrics import similarity as sim
perc_features_test.add_common_features()
perc_features_test.target_decoy_labels = perc_features_test.metadata["REVERSE"].apply(Percolator.get_target_decoy_label).to_numpy()

np.random.seed(1)
# add Prosit or Andromeda features
fragments_ratio = fr.FragmentsRatio(perc_features_test.pred_intensities, perc_features_test.true_intensities)
fragments_ratio.calc()
similarity = sim.SimilarityMetrics(perc_features_test.pred_intensities, perc_features_test.true_intensities, perc_features_test.mz)
similarity.calc(perc_features_test.all_features_flag)

perc_features_test.metrics_val = pd.concat(
    [perc_features_test.metrics_val, fragments_ratio.metrics_val, similarity.metrics_val], axis=1
)

lda_failed = False
perc_features_test.metrics_val[['spectral_angle', 'pearson_corr']]

idxs_below_lda_fdr = perc_features_test.apply_lda_and_get_indices_below_fdr(fdr_cutoff=perc_features_test.fdr_cutoff)
current_fdr = perc_features_test.fdr_cutoff
while len(idxs_below_lda_fdr) == 0:
    current_fdr += 0.01
    idxs_below_lda_fdr = perc_features_test.apply_lda_and_get_indices_below_fdr(fdr_cutoff=current_fdr)
    if current_fdr >= 0.1:
        lda_failed = True
        break

sampled_idxs = Percolator.sample_balanced_over_bins(perc_features_test.metadata[["RETENTION_TIME", "PREDICTED_IRT"]])

aligned_predicted_rts = Percolator.get_aligned_predicted_retention_times(
    perc_features_test.metadata["RETENTION_TIME"][sampled_idxs],
    perc_features_test.metadata["PREDICTED_IRT"][sampled_idxs],
    perc_features_test.metadata["PREDICTED_IRT"],
    "lowess",
)

perc_features_test.metrics_val["RT"] = perc_features_test.metadata["RETENTION_TIME"]
perc_features_test.metrics_val["pred_RT"] = perc_features_test.metadata["PREDICTED_IRT"]
perc_features_test.metrics_val["iRT"] = aligned_predicted_rts
perc_features_test.metrics_val["collision_energy_aligned"] = perc_features_test.metadata["COLLISION_ENERGY"] / 100.0
perc_features_test.metrics_val["abs_rt_diff"] = np.abs(perc_features_test.metadata["RETENTION_TIME"] - aligned_predicted_rts)

median_abs_error_lda_targets = np.median(perc_features_test.metrics_val["abs_rt_diff"].iloc[idxs_below_lda_fdr])
print(
    f"Median absolute error predicted vs observed retention time on targets < 1% FDR: {median_abs_error_lda_targets}"
)

perc_features_test.metadata["SCORE"]

perc_features_test.add_percolator_metadata_columns()

perc_features_test._reorder_columns_for_percolator()





import grpc
from tensorflow_serving.apis import get_model_metadata_pb2
from tensorflow_serving.apis import get_model_metadata_pb2_grpc

# Create a gRPC channel to the server
channel = grpc.insecure_channel('10.152.135.57:8500')

# Create a stub for the ModelService
stub = get_model_metadata_pb2_grpc.ModelServiceStub(channel)

# Create a GetModelMetadataRequest for all models
request = get_model_metadata_pb2.GetModelMetadataRequest()

# Send the request to the server and wait for the response
response = stub.GetModelMetadata(request)

# Parse the response to get the list of available models
model_names = []
for model_metadata in response.metadata:
    model_names.append(model_metadata.model_spec.name)

print("Available models:", model_names)

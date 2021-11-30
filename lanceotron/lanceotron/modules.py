from . import lanceotron as Ltron
from .utils import make_directory_name, build_model, calculate_pvalue_from_input
import os, sys
import numpy as np
import pyBigWig
import pickle
from sklearn.preprocessing import StandardScaler

# import tensorflow as tf
# from tensorflow import keras
import csv
import pkg_resources
import pandas as pd


def find_and_score_peaks(
    file: str,
    cutoff: float = 0,
    format: str = "web",
    threshold: int = 4,
    window: int = 400,
    folder: str = "./",
    skipheader: bool = False,
    **other,
) -> None:
    """Call peaks from a coverage track and score them using the LanceOTron model.

    Args:
        file (str): Path to bigwig track.
        threshold (int, optional): Initial threshold used for selecting candidate peaks. Defaults to 4.
        window (int, optional): Window size for rolling mean to use for selecting candidate peaks. Defaults to 400.
        folder (str, optional): Output folder. Defaults to "./".
        skipheader (bool, optional): Skip writing header. Defaults to False.
    """

    bigwig_file = file
    out_folder = make_directory_name(folder)
    initial_threshold = threshold

    min_peak_width = 50
    max_peak_width = 2000
    read_coverage_factor = 10 ** 9
    out_file_name = bigwig_file.split("/")[-1].split(".")[0] + "_L-tron.bed"

    pyBigWig_object = pyBigWig.open(bigwig_file)
    read_coverage_total = pyBigWig_object.header()["sumData"]
    read_coverage_rphm = read_coverage_total / read_coverage_factor
    pyBigWig_object.close()

    bigwig_data = Ltron.Bigwig_data(bigwig_file)
    genome_stats_dict = bigwig_data.get_genome_info()
    bed_file_out = []

    for chrom in genome_stats_dict:
        coverage_array_smooth = bigwig_data.make_chrom_coverage_map(
            genome_stats_dict[chrom], smoothing=window
        )
        enriched_region_coord_list = (
            Ltron.label_enriched_regions_dynamic_threshold_width(
                coverage_array_smooth,
                genome_stats_dict[chrom]["chrom_mean"] * initial_threshold,
                genome_stats_dict[chrom]["chrom_mean"],
                max_peak_width,
                min_region_size=min_peak_width,
            )
        )
        chrom_file_out = []
        if enriched_region_coord_list:
            wide_path = pkg_resources.resource_filename(
                "lanceotron.static", "standard_scaler_wide_v5_03.p"
            )
            deep_path = pkg_resources.resource_filename(
                "lanceotron.static", "standard_scaler_deep_v5_03.p"
            )

            coverage_array = (
                bigwig_data.make_chrom_coverage_map(genome_stats_dict[chrom])
                / read_coverage_rphm
            )
            X_wide_array, X_deep_array = Ltron.extract_signal_wide_and_deep_chrom(
                coverage_array, enriched_region_coord_list, read_coverage_rphm
            )
            standard_scaler_wide = pickle.load(open(wide_path, "rb"))
            X_wide_array_norm = standard_scaler_wide.transform(X_wide_array)
            X_wide_array_norm = np.expand_dims(X_wide_array_norm, axis=2)
            standard_scaler = StandardScaler()
            X_deep_array_norm_T = standard_scaler.fit_transform(X_deep_array.T)
            standard_scaler_deep = pickle.load(open(deep_path, "rb"))
            X_deep_array_norm = standard_scaler_deep.transform(X_deep_array_norm_T.T)
            X_deep_array_norm = np.expand_dims(X_deep_array_norm, axis=2)
            model = build_model()
            model_classifications = model.predict(
                [X_deep_array_norm, X_wide_array_norm], verbose=1
            )

            import tensorflow.keras.backend as K

            K.clear_session()

            for i, coord_pair in enumerate(enriched_region_coord_list):
                out_list = [
                    chrom,
                    coord_pair[0],
                    coord_pair[1],
                    model_classifications[0][i][0],
                    model_classifications[1][i][0],
                    model_classifications[2][i][0],
                ]
                X_wide_list = X_wide_array[i][:-1].tolist()
                X_wide_list = [100.0 if x > 10 else x for x in X_wide_list]
                out_list += X_wide_list
                chrom_file_out.append(out_list)
            bed_file_out += chrom_file_out

    output = pd.DataFrame(bed_file_out)
    output.columns = [
        "chrom",
        "start",
        "end",
        "overall_peak_score",
        "shape_score",
        "enrichment_score",
        "pvalue_chrom",
        "pvalue_10kb",
        "pvalue_20kb",
        "pvalue_30kb",
        "pvalue_40kb",
        "pvalue_50kb",
        "pvalue_60kb",
        "pvalue_70kb",
        "pvalue_80kb",
        "pvalue_90kb",
        "pvalue_100kb",
    ]
    output_thresholded = output.loc[output["overall_peak_score"] > cutoff, :]

    if format.lower() == "web":
        output_thresholded.to_csv(
            f"{out_folder}{out_file_name}", header=not skipheader, index=False, sep="\t"
        )
    elif format.lower() == "bed":
        output_thresholded.iloc[:, 0:3].to_csv(
            f"{out_folder}{out_file_name}", header=not skipheader, index=False, sep="\t"
        )
    else:
        raise NotImplementedError

    # with open(out_folder+out_file_name, 'w', newline='') as f:
    #    if not skipheader:
    #        f.write('chrom\tstart\tend\toverall_peak_score\tshape_score\tenrichment_score\tpvalue_chrom\tpvalue_10kb\tpvalue_20kb\tpvalue_30kb\tpvalue_40kb\tpvalue_50kb\tpvalue_60kb\tpvalue_70kb\tpvalue_80kb\tpvalue_90kb\tpvalue_100kb\n')
    #    bed_writer = csv.writer(f, delimiter='\t')
    #    import pdb; pdb.set_trace()
    #    bed_writer.writerows(bed_file_out)


def call_peaks_with_input(
    file: str,
    input: str,
    threshold: int = 4,
    window: int = 400,
    folder: str = "./",
    skipheader: bool = False,
    **other,
) -> None:
    """Call peaks from a coverage track with an input and score them using the LanceOTron model.

    Args:
        file (str): Path to bigwig track.
        input (str): Control input track used to calculate Poisson-based significance of peaks.
        threshold (int, optional): Initial threshold used for selecting candidate peaks. Defaults to 4.
        window (int, optional): Window size for rolling mean to use for selecting candidate peaks. Defaults to 400.
        folder (str, optional): Output folder. Defaults to "./".
        skipheader (bool, optional): Skip writing header. Defaults to False.
    """

    bigwig_file = file
    control_file = input
    out_folder = make_directory_name(folder)
    initial_threshold = threshold

    min_peak_width = 50
    max_peak_width = 2000
    read_coverage_factor = 10 ** 9
    out_file_name = bigwig_file.split("/")[-1].split(".")[0] + "_L-tron.bed"

    pyBigWig_object = pyBigWig.open(bigwig_file)
    read_coverage_total = pyBigWig_object.header()["sumData"]
    read_coverage_rphm = read_coverage_total / read_coverage_factor
    pyBigWig_object.close()

    bigwig_data = Ltron.Bigwig_data(bigwig_file)
    genome_stats_dict = bigwig_data.get_genome_info()
    bed_file_out = []

    for chrom in genome_stats_dict:
        coverage_array_smooth = bigwig_data.make_chrom_coverage_map(
            genome_stats_dict[chrom], smoothing=window
        )
        enriched_region_coord_list = (
            Ltron.label_enriched_regions_dynamic_threshold_width(
                coverage_array_smooth,
                genome_stats_dict[chrom]["chrom_mean"] * initial_threshold,
                genome_stats_dict[chrom]["chrom_mean"],
                max_peak_width,
                min_region_size=min_peak_width,
            )
        )
        chrom_file_out = []
        if enriched_region_coord_list:

            wide_path = pkg_resources.resource_filename(
                "lanceotron.static", "standard_scaler_wide_v5_03.p"
            )
            deep_path = pkg_resources.resource_filename(
                "lanceotron.static", "standard_scaler_deep_v5_03.p"
            )

            coverage_array = (
                bigwig_data.make_chrom_coverage_map(genome_stats_dict[chrom])
                / read_coverage_rphm
            )
            X_wide_array, X_deep_array = Ltron.extract_signal_wide_and_deep_chrom(
                coverage_array, enriched_region_coord_list, read_coverage_rphm
            )
            standard_scaler_wide = pickle.load(open(wide_path, "rb"))
            X_wide_array_norm = standard_scaler_wide.transform(X_wide_array)
            X_wide_array_norm = np.expand_dims(X_wide_array_norm, axis=2)
            standard_scaler = StandardScaler()
            X_deep_array_norm_T = standard_scaler.fit_transform(X_deep_array.T)
            standard_scaler_deep = pickle.load(open(deep_path, "rb"))
            X_deep_array_norm = standard_scaler_deep.transform(X_deep_array_norm_T.T)
            X_deep_array_norm = np.expand_dims(X_deep_array_norm, axis=2)
            model = build_model()
            model_classifications = model.predict(
                [X_deep_array_norm, X_wide_array_norm], verbose=1
            )
            pyBigWig_input = pyBigWig.open(control_file)
            read_coverage_total_input = pyBigWig_input.header()["sumData"]
            read_coverage_rphm_input = read_coverage_total_input / read_coverage_factor

            import tensorflow.keras.backend as K

            K.clear_session()
            for i, coord_pair in enumerate(enriched_region_coord_list):
                average_cov = (
                    coverage_array[coord_pair[0] : coord_pair[1]].mean()
                    * read_coverage_rphm
                )
                pvalue_input = calculate_pvalue_from_input(
                    chrom,
                    coord_pair[0],
                    coord_pair[1],
                    read_coverage_total,
                    read_coverage_total_input,
                    pyBigWig_input,
                    average_cov,
                )
                out_list = [
                    chrom,
                    coord_pair[0],
                    coord_pair[1],
                    model_classifications[0][i][0],
                    model_classifications[1][i][0],
                    model_classifications[2][i][0],
                    pvalue_input,
                ]
                X_wide_list = X_wide_array[i][:-1].tolist()
                X_wide_list = [100.0 if x > 10 else x for x in X_wide_list]
                out_list += X_wide_list
                chrom_file_out.append(out_list)
            pyBigWig_input.close()
            bed_file_out += chrom_file_out

    with open(out_folder + out_file_name, "w", newline="") as f:
        if not skipheader:
            f.write(
                "chrom\tstart\tend\toverall_peak_score\tshape_score\tenrichment_score\tpvalue_input\tpvalue_chrom\tpvalue_10kb\tpvalue_20kb\tpvalue_30kb\tpvalue_40kb\tpvalue_50kb\tpvalue_60kb\tpvalue_70kb\tpvalue_80kb\tpvalue_90kb\tpvalue_100kb\n"
            )
        bed_writer = csv.writer(f, delimiter="\t")
        bed_writer.writerows(bed_file_out)


def score_bed(
    file: str, bed: str, folder: str = "./", skipheader: bool = False, **other
) -> None:
    """Score an existing bed file using LanceOTron's model and a coverage track.

    Args:
        file (str): Input bigWig track.
        bed (str): Bed file of regions to score.
        folder (str): Output folder. Defaults to "./".
        skipheader (bool): Skip writing header. Defaults to False.
    """
    bigwig_file = file
    out_folder = make_directory_name(folder)
    bed_file = bed

    read_coverage_factor = 10 ** 9
    out_file_name = bigwig_file.split("/")[-1].split(".")[0] + "_L-tron.bed"

    pyBigWig_object = pyBigWig.open(bigwig_file)
    read_coverage_total = pyBigWig_object.header()["sumData"]
    read_coverage_rphm = read_coverage_total / read_coverage_factor
    pyBigWig_object.close()

    bed_list = Ltron.bed_file_to_list(bed_file)
    chroms_in_bed = []
    for bed_entry in bed_list:
        if bed_entry[0] not in chroms_in_bed:
            chroms_in_bed.append(bed_entry[0])

    bigwig_data = Ltron.Bigwig_data(bigwig_file)
    genome_stats_dict = bigwig_data.get_genome_info(include_special_chromosomes=True)
    bed_file_out = []

    for chrom in chroms_in_bed:
        enriched_region_coord_list = []
        for bed_entry in bed_list:
            if bed_entry[0] == chrom:
                enriched_region_coord_list.append([bed_entry[1], bed_entry[2]])
        chrom_file_out = []
        if enriched_region_coord_list:
            wide_path = pkg_resources.resource_filename(
                "lanceotron.static", "standard_scaler_wide_v5_03.p"
            )
            deep_path = pkg_resources.resource_filename(
                "lanceotron.static", "standard_scaler_deep_v5_03.p"
            )

            coverage_array = (
                bigwig_data.make_chrom_coverage_map(genome_stats_dict[chrom])
                / read_coverage_rphm
            )
            X_wide_array, X_deep_array = Ltron.extract_signal_wide_and_deep_chrom(
                coverage_array, enriched_region_coord_list, read_coverage_rphm
            )
            standard_scaler_wide = pickle.load(open(wide_path, "rb"))
            X_wide_array_norm = standard_scaler_wide.transform(X_wide_array)
            X_wide_array_norm = np.expand_dims(X_wide_array_norm, axis=2)
            standard_scaler = StandardScaler()
            X_deep_array_norm_T = standard_scaler.fit_transform(X_deep_array.T)
            standard_scaler_deep = pickle.load(open(deep_path, "rb"))
            X_deep_array_norm = standard_scaler_deep.transform(X_deep_array_norm_T.T)
            X_deep_array_norm = np.expand_dims(X_deep_array_norm, axis=2)
            model = build_model()
            model_classifications = model.predict(
                [X_deep_array_norm, X_wide_array_norm], verbose=1
            )

            import tensorflow.keras.backend as K

            K.clear_session()
            for i, coord_pair in enumerate(enriched_region_coord_list):
                out_list = [
                    chrom,
                    coord_pair[0],
                    coord_pair[1],
                    model_classifications[0][i][0],
                    model_classifications[1][i][0],
                    model_classifications[2][i][0],
                ]
                X_wide_list = X_wide_array[i][:-1].tolist()
                X_wide_list = [100.0 if x > 10 else x for x in X_wide_list]
                out_list += X_wide_list
                chrom_file_out.append(out_list)
            bed_file_out += chrom_file_out

    with open(out_folder + out_file_name, "w", newline="") as f:
        if not skipheader:
            f.write(
                "chrom\tstart\tend\toverall_peak_score\tshape_score\tenrichment_score\tpvalue_chrom\tpvalue_10kb\tpvalue_20kb\tpvalue_30kb\tpvalue_40kb\tpvalue_50kb\tpvalue_60kb\tpvalue_70kb\tpvalue_80kb\tpvalue_90kb\tpvalue_100kb\n"
            )
        bed_writer = csv.writer(f, delimiter="\t")
        bed_writer.writerows(bed_file_out)

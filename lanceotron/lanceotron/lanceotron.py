# lanceotron_classic_31

import numpy as np
import pyBigWig
from scipy.signal import savgol_filter
from scipy.stats import poisson
import csv


def average_array(array, window):
    if window > len(array):
        window = int(len(array) / 10)
        print(
            "averaging window larger than array, window adjusted to size {}".format(
                window
            )
        )
    if window == len(array):
        return_array = np.ones_like(array) * np.mean(array)
    else:
        return_array = np.cumsum(array, dtype=float)
        return_array[window:] = return_array[window:] - return_array[:-window]
        return_array /= window
        left_pad = int(window / 2)
        right_pad = window - left_pad
        return_array[left_pad:-right_pad] = return_array[window:]
        return_array[:left_pad] = np.mean(array[:left_pad])
        return_array[-right_pad:] = np.mean(array[-right_pad:])
    return return_array


def bed_file_to_list(bed_file, header=False):
    bed_list = []
    with open(bed_file, "r") as f:
        bed_reader = csv.reader(f, delimiter="\t")
        if header:
            next(bed_reader)
        for row in bed_reader:
            row_list = []
            for i, column in enumerate(row):
                if (i == 1) or (i == 2):
                    row_list.append(int(column))
                else:
                    row_list.append(column)
            bed_list.append(row_list)
    return bed_list


def calculate_max_local_mean(array, window_size_list):
    array_list = []
    for window in window_size_list:
        ave_array = average_array(array, window)
        array_list.append(ave_array)
    return np.maximum.reduce(array_list)


def calculate_min_local_mean(array, window_size_list):
    array_list = []
    for window in window_size_list:
        ave_array = average_array(array, window)
        array_list.append(ave_array)
    return np.minimum.reduce(array_list)


def label_enriched_regions_threshold(
    array, threshold, min_region_size=0, max_region_size=None
):
    truth_array = array > threshold
    enriched_region_list = []
    if truth_array.any:
        coord_array = np.where(truth_array[:-1] != truth_array[1:])[0]
        if coord_array.any:
            coord_array = coord_array + 1
            if len(coord_array) % 2 == 0:
                if truth_array[0] == True:
                    coord_array = np.insert(coord_array, 0, 0)
                    coord_array = np.append(coord_array, len(array))
            else:
                if truth_array[0] == True:
                    coord_array = np.insert(coord_array, 0, 0)
                else:
                    coord_array = np.append(coord_array, len(array))
            coord_array = coord_array.reshape(-1, 2)
            enriched_region_list = coord_array.tolist()
        else:
            if truth_array[0]:
                enriched_region_list = [[0, len(array)]]
    filtered_region_list = []
    if max_region_size is None:
        max_region_size = len(array)
    for coords in enriched_region_list:
        region_size = coords[1] - coords[0]
        if (region_size >= min_region_size) and (region_size <= max_region_size):
            filtered_region_list.append(coords)
    return filtered_region_list


def label_enriched_regions_dynamic_threshold_width(
    array,
    threshold,
    increment,
    max_width,
    min_region_size=0,
    return_intermediate_coords=False,
    recurse=False,
    level=0,
):
    level += 1
    enriched_region_list = label_enriched_regions_threshold(
        array, threshold, min_region_size
    )
    out_list = []
    new_regions = []
    corrected_list = []
    if enriched_region_list:
        for start_coord, end_coord in enriched_region_list:
            region_size = end_coord - start_coord
            if region_size > max_width:
                if level <= 50:
                    if return_intermediate_coords:
                        out_list.append([start_coord, end_coord])
                    new_regions = label_enriched_regions_dynamic_threshold_width(
                        array[start_coord:end_coord],
                        threshold + increment,
                        increment,
                        max_width,
                        min_region_size=min_region_size,
                        return_intermediate_coords=return_intermediate_coords,
                        recurse=True,
                        level=level,
                    )
                    if new_regions:
                        for region in new_regions:
                            corrected_list.append(
                                [region[0] + start_coord, region[1] + start_coord]
                            )
                    else:
                        out_list.append([start_coord, end_coord])
                else:
                    out_list.append([start_coord, end_coord])
            else:
                out_list.append([start_coord, end_coord])
    else:
        if recurse:
            out_list.append([0, len(array)])
    if new_regions:
        out_list.extend(corrected_list)
    return out_list


def select_best_savgol_window(array, window_list):
    lag_list = [x for x in range(1, 16)]
    window = 400
    fold_mean_enrichment = 2
    smooth_array = average_array(array, window)
    blank_signal_coord_list_filtered = []
    while not blank_signal_coord_list_filtered:
        enriched_regions_list = label_enriched_regions_threshold(
            smooth_array, fold_mean_enrichment * np.mean(array)
        )
        enriched_regions_array = np.array(enriched_regions_list).reshape(-1, 1)
        blank_signal_coord_list = []
        if len(enriched_regions_array) > 2:
            it = iter(enriched_regions_array[1:-1])
            for x in it:
                left_coord = x[0]
                right_coord = next(it)[0]
                region_size = right_coord - left_coord
                if (region_size > 2000) and (region_size < 10000):
                    blank_signal_coord_list.append(
                        [left_coord, right_coord, region_size]
                    )
        else:
            blank_signal_coord_list.append([0, len(array), len(array)])
        blank_signal_length = 0
        for coords in blank_signal_coord_list[:10000]:
            if np.count_nonzero(array[coords[0] : coords[1]]) > 100:
                blank_signal_coord_list_filtered.append(coords)
                blank_signal_length += coords[2]
        fold_mean_enrichment *= 2
    autocorr_list = []
    for lag in lag_list:
        autocorcoef_sum = 0
        for coords in blank_signal_coord_list_filtered:
            coeff = np.corrcoef(
                array[coords[0] : coords[1] - lag], array[coords[0] + lag : coords[1]]
            )[0][1]
            autocorcoef_sum += (coeff * coords[2]) / blank_signal_length
        autocorr_list.append(autocorcoef_sum)
    mean_autocorr = np.mean(autocorr_list)
    print(autocorr_list)
    print("mean coeffs:", mean_autocorr)
    autocorr_res_list = []
    for window_size in window_list:
        autocorcoef_sum = 0
        for coords in blank_signal_coord_list_filtered:
            savgol_window_test_array_smooth = savgol_filter(
                array[coords[0] : coords[1]], window_size, 2
            )
            savgol_window_test_array_res = (
                array[coords[0] : coords[1]] - savgol_window_test_array_smooth
            )
            input_padding = int((window_size - 1) / 2)
            coeff = np.corrcoef(
                savgol_window_test_array_res[input_padding : -input_padding - 1],
                savgol_window_test_array_res[input_padding + 1 : -input_padding],
            )[0][1]
            autocorcoef_sum += (coeff * coords[2]) / blank_signal_length
        autocorr_res_list.append(autocorcoef_sum)
        print(window_size, ":", autocorcoef_sum)
    window_list_index = (np.abs(autocorr_res_list - mean_autocorr)).argmin()
    print(window_list[window_list_index])
    return window_list[window_list_index]


def label_enriched_regions_savgol(
    array,
    window_len,
    min_region_size,
    polynomial=2,
    window_list=None,
    max_region_size=None,
    threshold=0,
):
    if window_len == "auto":
        if window_list:
            window_len = select_best_savgol_window(array, window_list)
        else:
            window_list = [51, 101, 201, 401, 801]
            window_len = select_best_savgol_window(array, window_list)
            print(window_len)
    second_deriv_array = savgol_filter(array, window_len, 2, deriv=2)
    zero_crossings = np.where(np.diff(np.signbit(second_deriv_array)))[0]
    enriched_region_list = []
    if zero_crossings.any:
        if len(zero_crossings) % 2 == 0:
            if second_deriv_array[zero_crossings[0]] < 0:
                zero_crossings = np.insert(zero_crossings, 0, 0)
                zero_crossings = np.append(zero_crossings, len(second_deriv_array))
        else:
            if second_deriv_array[zero_crossings[0]] > 0:
                zero_crossings = np.append(zero_crossings, len(second_deriv_array))
            else:
                zero_crossings = np.insert(zero_crossings, 0, 0)
        zero_crossings = zero_crossings.reshape(-1, 2)
        enriched_region_list = zero_crossings.tolist()
    filtered_region_list = []
    if max_region_size is None:
        max_region_size = len(array)
    for coords in enriched_region_list:
        region_size = coords[1] - coords[0]
        max_height = np.amax(array[coords[0] : coords[1]])
        if (
            (region_size >= min_region_size)
            and (region_size <= max_region_size)
            and (max_height > threshold)
        ):
            filtered_region_list.append(coords)
    return filtered_region_list


def extract_signal_single(
    pyBigWig_object,
    chrom_name,
    start_coord,
    end_coord,
    chrom_len,
    features,
    zero_pad=True,
):
    X = np.zeros(features)
    region_length = end_coord - start_coord
    if ((start_coord + int(region_length / 2) - int(features / 2)) >= 0) and (
        (start_coord + int(region_length / 2) - int(features / 2) + features)
        <= chrom_len
    ):
        if zero_pad:
            if region_length > features:
                region_start = start_coord + int(region_length / 2) - int(features / 2)
                X = np.nan_to_num(
                    np.array(
                        pyBigWig_object.values(
                            chrom_name, region_start, region_start + features
                        )
                    )
                )
            else:
                pad_length = int((features - region_length) / 2)
                X[pad_length : pad_length + region_length] = np.nan_to_num(
                    np.array(pyBigWig_object.values(chrom_name, start_coord, end_coord))
                )
        else:
            region_start = start_coord + int(region_length / 2) - int(features / 2)
            X = np.nan_to_num(
                np.array(
                    pyBigWig_object.values(
                        chrom_name, region_start, region_start + features
                    )
                )
            )
    return X


def calculate_poisson_pvalue_array(enrichment_list, lambda_list):
    pvalue_array = np.zeros((len(enrichment_list), len(lambda_list)))
    for i, enrichment in enumerate(enrichment_list):
        for j, lambda_val in enumerate(lambda_list):
            with np.errstate(divide="ignore"):
                pvalue_array[i][j] = -1 * np.log10(
                    1 - poisson.cdf(enrichment, lambda_val)
                )
    pvalue_array = np.reshape(pvalue_array, (1, -1))
    pvalue_array_infs_indicies = np.where(np.isinf(pvalue_array))
    pvalue_array[pvalue_array_infs_indicies] = 999
    return pvalue_array


def extract_signal_wide_and_deep_single_entry(
    pyBigWig_object,
    chrom_name,
    start_coord,
    end_coord,
    chrom_len,
    chrom_mean,
    seq_depth,
):
    features = 2000
    X_deep = np.zeros(features)
    lambda_cov = 100000
    lambda_count = 10
    pvalue_array = np.zeros(lambda_count + 1)
    region_length = end_coord - start_coord
    interval = int((lambda_cov / lambda_count))
    length_list = [interval * (x + 1) for x in range(lambda_count)]
    max_height = (
        pyBigWig_object.stats(
            chrom_name, start_coord, end_coord, type="max", exact=True
        )[0]
        / seq_depth
    )
    if max_height > 0:
        region_start = start_coord + int(region_length / 2) - int(features / 2)
        if (region_start >= 0) and (region_start + features <= chrom_len):
            X_deep = (
                np.nan_to_num(
                    np.array(
                        pyBigWig_object.values(
                            chrom_name, region_start, region_start + features
                        )
                    )
                )
            ) / seq_depth
        # logic for edge cases - if not edge, pull out whole area and use array slices to find means
        # else pull out larger and larger pieces:
        #     if beyond bound two-sided use last mean; one sided, use signal until bound
        lambda_list = [chrom_mean / seq_depth]
        if (start_coord + int(region_length / 2) - lambda_cov >= 0) and (
            start_coord + int(region_length / 2) + lambda_cov <= chrom_len
        ):
            lambda_cov_array = (
                np.nan_to_num(
                    np.array(
                        pyBigWig_object.values(
                            chrom_name,
                            start_coord + int(region_length / 2) - lambda_cov,
                            start_coord + int(region_length / 2) + lambda_cov,
                        )
                    )
                )
                / seq_depth
            )
            for length in length_list:
                pad = lambda_cov - length
                if pad == 0:
                    lambda_list.append(np.mean(lambda_cov_array))
                else:
                    lambda_list.append(np.mean(lambda_cov_array[pad:-pad]))
        else:
            for i, length in enumerate(length_list):
                if (start_coord + int(region_length / 2) - length >= 0) or (
                    start_coord + int(region_length / 2) + length <= chrom_len
                ):
                    if start_coord + int(region_length / 2) - length >= 0:
                        start_coord_for_mean = (
                            start_coord + int(region_length / 2) - length
                        )
                    else:
                        start_coord_for_mean = 0
                    if start_coord + int(region_length / 2) + length <= chrom_len:
                        end_coord_for_mean = start_coord + int(region_length / 2)
                    else:
                        end_coord_for_mean = chrom_len
                    lambda_entry = pyBigWig_object.stats(
                        chrom_name,
                        start_coord_for_mean,
                        end_coord_for_mean,
                        type="mean",
                        exact=True,
                    )
                    lambda_list.append(lambda_entry)
                else:
                    lambda_list.append(lambda_list[i])
        enrichment_list = [max_height]
        pvalue_array = calculate_poisson_pvalue_array(enrichment_list, lambda_list)
    X_wide = np.append(pvalue_array, seq_depth)
    return X_wide, X_deep


def extract_signal_wide_and_deep_chrom(chrom_cov_array, coord_list, seq_depth):
    features = 2000
    X_deep = np.zeros((len(coord_list), features))
    lambda_cov = 100000
    lambda_count = 10
    interval = int((lambda_cov / lambda_count))
    length_list = [interval * (x + 1) for x in range(lambda_count)]
    X_wide = np.zeros((len(coord_list), lambda_count + 2))
    chrom_len = len(chrom_cov_array)
    chrom_mean = chrom_cov_array.mean()
    for i, coord in enumerate(coord_list):
        pvalue_array = np.zeros(lambda_count + 1)
        start_coord = coord[0]
        end_coord = coord[1]
        region_length = end_coord - start_coord
        max_height = chrom_cov_array[start_coord:end_coord].max()
        if max_height > 0:
            region_start = start_coord + int(region_length / 2) - int(features / 2)
            coverage_length = len(
                chrom_cov_array[region_start : region_start + features]
            )
            if coverage_length == features:
                X_deep[i] = chrom_cov_array[region_start : region_start + features]
            else:
                if chrom_len < features:
                    pad = int((features - chrom_len) / 2)
                    X_deep[i][pad : pad + chrom_len] = chrom_cov_array
                elif region_start < 0:
                    region_start = 0
                    X_deep[i] = chrom_cov_array[region_start : region_start + features]
                else:
                    region_start = chrom_len - features
                    X_deep[i] = chrom_cov_array[region_start : region_start + features]
            # logic for edge cases - if not edge, pull out whole area and use array slices to find means
            # else pull out larger and larger pieces:
            #     if beyond bound two-sided use last mean; one sided, use signal until bound
            lambda_list = [chrom_mean / seq_depth]
            if (start_coord + int(region_length / 2) - lambda_cov >= 0) and (
                start_coord + int(region_length / 2) + lambda_cov <= chrom_len
            ):
                lambda_cov_array = chrom_cov_array[
                    start_coord
                    + int(region_length / 2)
                    - lambda_cov : start_coord
                    + int(region_length / 2)
                    + lambda_cov
                ]
                for length in length_list:
                    pad = lambda_cov - length
                    if pad == 0:
                        lambda_list.append(lambda_cov_array.mean())
                    else:
                        lambda_list.append(lambda_cov_array[pad:-pad].mean())
            else:
                for j, length in enumerate(length_list):
                    if (start_coord + int(region_length / 2) - length >= 0) or (
                        start_coord + int(region_length / 2) + length <= chrom_len
                    ):
                        if start_coord + int(region_length / 2) - length >= 0:
                            start_coord_for_mean = (
                                start_coord + int(region_length / 2) - length
                            )
                        else:
                            start_coord_for_mean = 0
                        if start_coord + int(region_length / 2) + length <= chrom_len:
                            end_coord_for_mean = start_coord + int(region_length / 2)
                        else:
                            end_coord_for_mean = chrom_len
                        lambda_entry = chrom_cov_array[
                            start_coord_for_mean:end_coord_for_mean
                        ].mean()
                        lambda_list.append(lambda_entry)
                    else:
                        lambda_list.append(lambda_list[j])
            enrichment_list = [max_height]
            pvalue_array = calculate_poisson_pvalue_array(enrichment_list, lambda_list)
        X_wide[i] = np.append(pvalue_array, seq_depth)
    return X_wide, X_deep


def extract_signal_single_dl_input1(
    pyBigWig_object, chrom_name, start_coord, end_coord, features
):
    X = np.zeros((features, 3))
    region_length = end_coord - start_coord
    if region_length > features:
        region_start = start_coord + int(region_length / 2) - int(features / 2)
        X[:, 0] = np.nan_to_num(
            np.array(
                pyBigWig_object.values(
                    chrom_name, region_start, region_start + features
                )
            )
        )
        X[:, 1] = np.nan_to_num(
            np.array(
                pyBigWig_object.values(
                    chrom_name, region_start - features, region_start
                )
            )
        )
        X[:, 2] = np.nan_to_num(
            np.array(
                pyBigWig_object.values(
                    chrom_name, region_start + features, region_start + (2 * features)
                )
            )
        )
    else:
        pad_length = int((features - region_length) / 2)
        X[pad_length : pad_length + region_length, 0] = np.nan_to_num(
            np.array(pyBigWig_object.values(chrom_name, start_coord, end_coord))
        )
        X[:, 1] = np.nan_to_num(
            np.array(
                pyBigWig_object.values(chrom_name, start_coord - features, start_coord)
            )
        )
        X[:, 2] = np.nan_to_num(
            np.array(
                pyBigWig_object.values(chrom_name, end_coord, end_coord + features)
            )
        )
    return X


def extract_signal_single_dl_input2_old(
    pyBigWig_object,
    chrom_name,
    start_coord,
    end_coord,
    chrom_mean,
    chrom_std,
    kb_to_collect_stats=20,
):
    X = np.zeros((kb_to_collect_stats, 2))
    mid = int((end_coord - start_coord) / 2) + start_coord
    full_range = int(kb_to_collect_stats * 1000)
    half_range = int(full_range / 2)
    region_coverage = np.nan_to_num(
        np.array(
            pyBigWig_object.values(
                chrom_name, (mid - half_range), (mid - half_range) + full_range
            )
        )
    )
    region_coverage = (region_coverage - chrom_mean) / chrom_std

    X[0, 0] = np.mean(
        region_coverage[
            start_coord - (mid - half_range) : end_coord - (mid - half_range)
        ]
    )
    X[0, 1] = np.std(
        region_coverage[
            start_coord - (mid - half_range) : end_coord - (mid - half_range)
        ]
    )

    for i, region_length in enumerate(
        range(1000, (1000 * kb_to_collect_stats) + 1000, 1000)
    ):
        if i > 0:
            half_length = int(region_length / 2)
            X[i, 0] = np.mean(
                region_coverage[
                    half_range
                    - half_length : (half_range - half_length)
                    + region_length
                ]
            )
            X[i, 1] = np.std(
                region_coverage[
                    half_range
                    - half_length : (half_range - half_length)
                    + region_length
                ]
            )
    return X


def extract_signal_single_dl_input2(
    pyBigWig_object,
    chrom_name,
    start_coord,
    end_coord,
    chrom_len,
    chrom_mean,
    chrom_std,
    seq_depth,
    kb_to_collect_stats=20,
):
    X = np.zeros((kb_to_collect_stats + 2, 2))
    mid = int((end_coord - start_coord) / 2) + start_coord
    full_range = int(kb_to_collect_stats * 1000)
    half_range = int(full_range / 2)
    if ((mid - half_range) >= 0) and (((mid - half_range) + full_range) <= chrom_len):
        region_coverage = np.nan_to_num(
            np.array(
                pyBigWig_object.values(
                    chrom_name, (mid - half_range), (mid - half_range) + full_range
                )
            )
        )
    else:
        region_coverage = np.zeros(full_range)
        # region_coverage[start_coord-(mid-half_range):end_coord-(mid-half_range)] = np.nan_to_num(np.array(pyBigWig_object.values(chrom_name, start_coord, end_coord)))
    X[0, 0] = np.mean(
        region_coverage[
            start_coord - (mid - half_range) : end_coord - (mid - half_range)
        ]
    )
    X[0, 1] = np.std(
        region_coverage[
            start_coord - (mid - half_range) : end_coord - (mid - half_range)
        ]
    )
    X[1, 0] = chrom_mean
    X[1, 1] = chrom_std

    for i, region_length in enumerate(
        range(1000, (1000 * kb_to_collect_stats) + 1000, 1000)
    ):
        half_length = int(region_length / 2)
        X[i + 2, 0] = np.mean(region_coverage[region_length - 1000 : region_length])
        X[i + 2, 1] = np.std(region_coverage[region_length - 1000 : region_length])
    X = np.reshape(X, (1, -1))
    X = np.append(X, seq_depth)
    return X


def find_summits(coverage_array, enriched_region_coord_list):
    summit_list = []
    for region in enriched_region_coord_list:
        region_array = coverage_array[region[0]:region[1]]
        max_index = np.argmax(region_array)
        summit_list.append(region[0] + max_index)
    return summit_list


class Bigwig_data:
    import pyBigWig
    import numpy as np
    from scipy.signal import savgol_filter

    def __init__(self, bigWig_file):
        self.bigWig_file = bigWig_file

    def get_chrom_info(self, chrom_name):
        pyBigWig_object = pyBigWig.open(self.bigWig_file)
        chrom_stats_dict = {
            "chrom_name": chrom_name,
            "chrom_len": pyBigWig_object.chroms(chrom_name),
            "chrom_mean": pyBigWig_object.stats(chrom_name, type="mean", exact=True)[0],
            "chrom_std": pyBigWig_object.stats(chrom_name, type="std", exact=True)[0],
        }
        pyBigWig_object.close()
        return chrom_stats_dict

    def get_genome_info(self, include_special_chromosomes=False):
        genome_stats_dict = {}
        chrom_list = []
        pyBigWig_object = pyBigWig.open(self.bigWig_file)
        for chrom_name in pyBigWig_object.chroms():
            if include_special_chromosomes:
                chrom_list.append(chrom_name)
            else:
                if (
                    (not chrom_name.startswith("chrUn"))
                    and ("_" not in chrom_name)
                    and (chrom_name != "chrM")
                    and (chrom_name != "chrEBV")
                ):
                    chrom_list.append(chrom_name)
        pyBigWig_object.close()
        for chrom_name in chrom_list:
            chrom_stats_dict = self.get_chrom_info(chrom_name)
            genome_stats_dict[chrom_name] = chrom_stats_dict
        return genome_stats_dict

    def make_chrom_coverage_map(self, chrom_stats_dict, smoothing=None, savgol=False):
        chrom_coverage_array = np.zeros(chrom_stats_dict["chrom_len"])
        pyBigWig_object = pyBigWig.open(self.bigWig_file)
        if pyBigWig_object.intervals(chrom_stats_dict["chrom_name"]) is not None:
            for reads in pyBigWig_object.intervals(chrom_stats_dict["chrom_name"]):
                chrom_coverage_array[reads[0] : reads[1]] = reads[2]
        else:
            print("no reads found on this chromosome")
        pyBigWig_object.close()
        if smoothing is None:
            return chrom_coverage_array
        else:
            if savgol:
                smooth_array = savgol_filter(chrom_coverage_array, smoothing, 2)
                return smooth_array
            else:
                smooth_array = average_array(chrom_coverage_array, smoothing)
                return smooth_array

    def get_chrom_info_make_coverage_map(
        self, chrom_name, smoothing=None, savgol=False, return_chrom_stats_dict=False
    ):
        chrom_stats_dict = self.get_chrom_info(chrom_name)
        chrom_coverage_array = self.make_chrom_coverage_map(
            chrom_stats_dict, smoothing, savgol
        )
        if return_chrom_stats_dict:
            return chrom_coverage_array, chrom_stats_dict
        else:
            return chrom_coverage_array

    def make_chrom_coverage_map_label_enriched_regions_threshold(
        self,
        chrom_stats_dict,
        threshold,
        min_region_size=0,
        max_region_size=None,
        smoothing=None,
        savgol=False,
        return_chrom_coverage_array=False,
    ):
        chrom_coverage_array = self.make_chrom_coverage_map(
            chrom_stats_dict, smoothing, savgol
        )
        enriched_regions_list = label_enriched_regions_threshold(
            chrom_coverage_array, threshold, min_region_size, max_region_size
        )
        if return_chrom_coverage_array:
            return enriched_regions_list, chrom_coverage_array
        else:
            return enriched_regions_list

    def make_chrom_coverage_map_label_enriched_regions_savgol(
        self,
        chrom_stats_dict,
        smoothing="auto",
        min_region_size=0,
        max_region_size=None,
        threshold=0,
        return_chrom_coverage_array=False,
    ):
        chrom_coverage_array = self.make_chrom_coverage_map(chrom_stats_dict)
        enriched_regions_list = label_enriched_regions_savgol(
            chrom_coverage_array,
            smoothing,
            min_region_size,
            max_region_size,
            threshold=threshold,
        )
        if return_chrom_coverage_array:
            return enriched_regions_list, chrom_coverage_array
        else:
            return enriched_regions_list

    def get_chrom_info_make_chrom_coverage_map_label_enriched_regions_threshold(
        self,
        chrom_name,
        threshold,
        min_region_size=0,
        max_region_size=None,
        smoothing=None,
        savgol=False,
        return_chrom_coverage_array_stats_dict=False,
    ):
        if return_chrom_coverage_array_stats_dict:
            (
                chrom_coverage_array,
                chrom_stats_dict,
            ) = self.get_chrom_info_make_coverage_map(
                chrom_name, smoothing, savgol, return_chrom_stats_dict=True
            )
        else:
            chrom_coverage_array = self.get_chrom_info_make_coverage_map(
                chrom_name, smoothing, savgol
            )
        enriched_regions_list = label_enriched_regions_threshold(
            chrom_coverage_array, threshold, min_region_size, max_region_size
        )
        if return_chrom_coverage_array_stats_dict:
            return enriched_regions_list, chrom_coverage_array, chrom_stats_dict
        else:
            return enriched_regions_list

    def get_chrom_info_make_chrom_coverage_map_label_enriched_regions_savgol(
        self,
        chrom_name,
        smoothing="auto",
        min_region_size=0,
        max_region_size=None,
        threshold=0,
        return_chrom_coverage_array=False,
    ):
        chrom_coverage_array = self.get_chrom_info_make_coverage_map(chrom_name)
        enriched_regions_list = label_enriched_regions_savgol(
            chrom_coverage_array, smoothing, min_region_size, max_region_size
        )
        if return_chrom_coverage_array:
            return enriched_regions_list, chrom_coverage_array
        else:
            return enriched_regions_list

    def extract_signal_chrom(
        self, chrom_stats_dict, coord_list, features, normalize=True, zero_pad=True
    ):
        X = np.zeros((len(coord_list), features))
        pyBigWig_object = pyBigWig.open(self.bigWig_file)
        for i, coord in enumerate(coord_list):
            region_length = coord[1] - coord[0]
            if zero_pad:
                if region_length > features:
                    region_start = coord[0] + int(region_length / 2) - int(features / 2)
                    X[i] = np.nan_to_num(
                        np.array(
                            pyBigWig_object.values(
                                chrom_stats_dict["chrom_name"],
                                region_start,
                                region_start + features,
                            )
                        )
                    )
                else:
                    pad_length = int((features - region_length) / 2)
                    X[i, pad_length : pad_length + region_length] = np.nan_to_num(
                        np.array(
                            pyBigWig_object.values(
                                chrom_stats_dict["chrom_name"], coord[0], coord[1]
                            )
                        )
                    )
            else:
                region_start = coord[0] + int(region_length / 2) - int(features / 2)
                X[i] = np.nan_to_num(
                    np.array(
                        pyBigWig_object.values(
                            chrom_stats_dict["chrom_name"],
                            region_start,
                            region_start + features,
                        )
                    )
                )
        if normalize:
            X = (X - chrom_stats_dict["chrom_mean"]) / chrom_stats_dict["chrom_std"]
        return X

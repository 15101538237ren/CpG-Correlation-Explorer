# -*- coding: utf-8 -*-
import os, pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr

data_dir = "../data"
chromHMM_labels = range(1, 16)#
distance_range = range(2, 1000 + 1)
chromosomes = ["chr%d" % i for i in range(1, 2)]
FEATURE_NAME = "feature"
CORR_NAME = "corr"
def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def step1_split_bed_file_by_chromHMM_label(input_bed_fp, outdir, col_index=-1):
    mkdir(outdir)
    try:
        fmt = "%s\t%d\t%d" + "".join(["\t%.4f" for item in range(15)]) + "\t%d"
        df = pd.read_csv(input_bed_fp, sep='\t', header=None).values
        for chrHMMlbl in chromHMM_labels:
            df_lines = df[df[:, col_index] == chrHMMlbl, ]
            output_fp = os.path.join(outdir, "%s.bed" % chrHMMlbl)
            np.savetxt(output_fp, df_lines[:], delimiter='\n', fmt=fmt)
            print("save %s sucessful" % chrHMMlbl)
    except Exception as e:
        print(e)

def partition_CpG_pairs_dict_into_training_validation_and_testing_set(CpG_pair_withd_chrHMM, training_ratio):
    CpG_pair_withd_chrHMM = CpG_pair_withd_chrHMM[0]

    training_set = {chrHMMlbl: {d: {FEATURE_NAME:[]} for d in distance_range} for chrHMMlbl in chromHMM_labels}
    testing_set = {chrHMMlbl: {d: {FEATURE_NAME:[]} for d in distance_range} for chrHMMlbl in chromHMM_labels}

    for chrHMMlbl in chromHMM_labels:
        for d in distance_range:
            feature_vec_list = CpG_pair_withd_chrHMM[chrHMMlbl][d][FEATURE_NAME]
            if len(feature_vec_list):
                training_set[chrHMMlbl][d][FEATURE_NAME], testing_set[chrHMMlbl][d][FEATURE_NAME] = train_test_split(feature_vec_list, train_size=training_ratio)
    return [training_set], [testing_set]

def step2_construct_cpg_pairs_with_d_distance(input_dir, load = False):
    pkl_fp = os.path.join(data_dir, "CpG_pairs_with_d_ChromHMM.pkl")
    if load:
        with open(pkl_fp, 'rb') as pkl:
            # dump information to that file
            CpG_pair_withd_chrHMM = pickle.load(pkl)
    else:
        CpG_pair_withd_chrHMM = {} # the dict for storing the data of CpG pairs with d distance, for each chromHMM type.
        for chrHMMlbl in chromHMM_labels:
            CpG_pair_withd_chrHMM[chrHMMlbl] = {d: {FEATURE_NAME:[]} for d in distance_range}
            input_bed_fp = os.path.join(input_dir, "%s.bed" % chrHMMlbl)
            df = pd.read_csv(input_bed_fp, sep='\t', header=None).values
            for chrm in chromosomes:
                df_lines = df[df[:, 0] == chrm, :]
                CpG_pos_and_data_dict = {}
                for item in df_lines:
                    pos = item[1]
                    feature_vec = item[3:]
                    CpG_pos_and_data_dict[pos] = feature_vec
                for d in distance_range:
                    print("Step 2: HMM %d, %s, d = %d" % (chrHMMlbl, chrm, d))
                    for pos in CpG_pos_and_data_dict.keys():
                        pos_off_by_d = pos + d
                        if pos_off_by_d in CpG_pos_and_data_dict.keys():
                            kdiff = abs(CpG_pos_and_data_dict[pos][1] - CpG_pos_and_data_dict[pos_off_by_d][1])
                            feature_of_CpG_pairs = [CpG_pos_and_data_dict[pos], CpG_pos_and_data_dict[pos_off_by_d], kdiff]
                            CpG_pair_withd_chrHMM[chrHMMlbl][d][FEATURE_NAME].append(feature_of_CpG_pairs)
        #training_set, testing_set = partition_CpG_pairs_dict_into_training_validation_and_testing_set([CpG_pair_withd_chrHMM], training_ratio=0.5)
        CpG_pair_withd_chrHMM = [CpG_pair_withd_chrHMM]
        with open(pkl_fp, 'wb') as pkl:
            # dump information to that file
            pickle.dump(CpG_pair_withd_chrHMM, pkl, -1)
    return CpG_pair_withd_chrHMM

def local_feature_keys(d, local_radius):
    if (d - local_radius) not in distance_range:
        ds = range(d, d + 2 * local_radius + 1)
    elif (d + local_radius + 1) not in distance_range:
        ds = range(d - 2 * local_radius + 1, d)
    else:
        ds = range(d - local_radius, d + local_radius + 1)
    return ds

def calc_local_correlation(dataset, local_radius, target_feature_col = 5, NUM_MIN_SAMPLE = 5):
    dataset = dataset[0]
    COL_INDEX_IN_DF = target_feature_col - 3 - 1
    for chrHMMlbl in chromHMM_labels:
        for d in distance_range:
            print("Step 3: HMM %d, d = %d" % (chrHMMlbl, d))
            ds = local_feature_keys(d, local_radius)
            vec_list = []
            for dd in ds:
                feature_data = dataset[chrHMMlbl][dd][FEATURE_NAME]
                if feature_data:
                    for feature in feature_data:
                        vec_list.append([feature[0][COL_INDEX_IN_DF], feature[1][COL_INDEX_IN_DF]])
            if len(vec_list) > NUM_MIN_SAMPLE:
                vec_list = np.array(vec_list)
                corr, _ = pearsonr(vec_list[:, 0], vec_list[:, 1])
                dataset[chrHMMlbl][d][CORR_NAME] = corr
            else:
                dataset[chrHMMlbl][d][CORR_NAME] = -2
    return [dataset]

def merge_diff_chromHMM(dataset):
    dataset = dataset[0]
    new_dataset = {d: [] for d in distance_range}
    for d in distance_range:
        for chrHMMlbl in chromHMM_labels:
            corr = dataset[chrHMMlbl][d][CORR_NAME]
            if corr > -2:
                feature_list = dataset[chrHMMlbl][d][FEATURE_NAME]
                for feature in feature_list:
                    new_dataset[d].append([feature, corr])
    return new_dataset
def step3_get_local_correlation_for_each_CpG_pair(local_radius= 2, load = False, training_set= None, testing_set = None):
    pkl_fp = os.path.join(data_dir, "training_testing_corr_before_merging.pkl")
    if load:
        with open(pkl_fp, 'rb') as pkl:
            [training_set, testing_set] = pickle.load(pkl)
    else:
        training_set = calc_local_correlation(training_set, local_radius)
        testing_set = calc_local_correlation(testing_set, local_radius)

        with open(pkl_fp, 'wb') as pkl:
            pickle.dump([training_set, testing_set], pkl, -1)
    # training_set = merge_diff_chromHMM(training_set)
    # testing_set = merge_diff_chromHMM(testing_set)
    return training_set, testing_set

if __name__ == "__main__":
    input_bed_fp = os.path.join(data_dir, "K-100bp-CpGd_HMM.bed")
    outdir = os.path.join(data_dir, "K-100bp-CpGd_split")
    #step1_split_bed_file_by_chromHMM_label(input_bed_fp, outdir)

    CpG_pair_withd_chrHMM = step2_construct_cpg_pairs_with_d_distance(outdir, load= False)
    #training_set, testing_set = step3_get_local_correlation_for_each_CpG_pair(load=True)
#
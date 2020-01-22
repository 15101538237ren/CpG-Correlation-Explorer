# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_dir = "../data"
chromHMM_labels = range(1, 16)
distance_range = range(2, 1000 + 1)
chromosomes = ["chr%d" % i for i in range(1, 23)]

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def step1_split_bed_file_by_chromHMM_label(input_bed_fp, outdir, col_index=-1):
    mkdir(outdir)
    try:
        fmt = "%s\t%d\t%d" + "".join(["\t%.4f" for item in range(15)])
        df = pd.read_csv(input_bed_fp, sep='\t', header=None).values
        for chrHMMlbl in chromHMM_labels:
            df_lines = df[df[:, col_index] == chrHMMlbl, 0: col_index]
            output_fp = os.path.join(outdir, "%s.bed" % chrHMMlbl)
            np.savetxt(output_fp, df_lines[:], delimiter='\n', fmt=fmt)
            print("save %s sucessful" % chrHMMlbl)
    except Exception as e:
        print(e)
def step2_construct_cpg_pairs_with_d_distance(input_dir):
    CpG_pair_withd_chrHMM = {} # the dict for storing the data of CpG pairs with d distance, for each chromHMM type.
    for chrHMMlbl in chromHMM_labels:
        CpG_pair_withd_chrHMM[chrHMMlbl] = {d: [] for d in distance_range}
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
                for pos in CpG_pos_and_data_dict.keys():
                    pos_off_by_d = pos + d
                    if pos_off_by_d in CpG_pos_and_data_dict.keys():
                        feature_of_CpG_pairs = [CpG_pos_and_data_dict[pos], CpG_pos_and_data_dict[pos_off_by_d]]
                        CpG_pair_withd_chrHMM[chrHMMlbl][d].append(feature_of_CpG_pairs)

if __name__ == "__main__":

    input_bed_fp = os.path.join(data_dir, "K-100bp-CpGd_HMM.bed")
    outdir = os.path.join(data_dir, "K-100bp-CpGd_split")
    #step1_split_bed_file_by_chromHMM_label(input_bed_fp, outdir)

    step2_construct_cpg_pairs_with_d_distance(outdir)
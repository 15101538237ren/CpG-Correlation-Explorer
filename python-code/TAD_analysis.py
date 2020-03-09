# -*- coding: utf-8 -*-
import os, random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from localreg import *
import seaborn as sns
from bisect import *

sns.set(color_codes=True)
data_dir = "../data"
figure_dir = "../figures/corr_around_TADs"
EACH_SUB_FIG_SIZE = 5
nchr = 22
window_sz = 2000
min_nCpG = 5  # number of the minimal CpGs number for the regions in TAD boundaries.


def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def read_TAD_k_intersection_file(input_fp, isTAD=False):
    df = pd.read_csv(input_fp, sep="\t", header=None).values


def histogram_of_TAD_length(TAD_fp, non_TAD_fp, fig_format="png"):
    NBIN = 30
    cm = plt.get_cmap('gist_rainbow')
    mkdir(figure_dir)
    fig_fp = os.path.join(figure_dir, "TAD_%dbp_updownstream.%s" % (window_sz,fig_format))
    N_COL = 2
    N_ROW = 2
    fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))

    ax = axs[0][0]  # Histogram of TAD size
    tad_fp = os.path.join(data_dir, "TADs", "h1_all.hg19.20kby40k.all.finaldomaincalls.bed")
    df = pd.read_csv(tad_fp, sep="\t", header=None).values
    tad_sizes = np.abs(np.array(df[:, -1] - df[:, -2]).astype(float))
    _ = ax.hist(tad_sizes, bins=NBIN, density=True, color='blue', edgecolor='black',
                alpha=0.5, linewidth=0.5)
    ax.set_xlim([0, 5e6])
    ax.set_ylabel('Probability')
    ax.set_xlabel('Genomic Distance(bp)')
    ax.set_title("Histogram of TAD sizes of hESC")

    [x, y1, y2] = correlation_between_k_in_TAD_bounaries_and_k_in_non_TAD_large_slide_wd(TAD_fp, non_TAD_fp)
    ax = axs[0][1]  # Histogram of TAD size
    ax.scatter(x, y1, s=2, color="red")
    ax.scatter(x, y2, s=2, color="blue")
    ax.set_ylim([0, 1])
    ax.set_ylabel('|Correlation|')
    ax.set_xlabel('Genomic Distance(bp)')
    ax.set_title("Correlation k in TADs vs not")

    ax = axs[1][0]  # Histogram of TAD size
    my1 = np.mean(y1[~np.isnan(y1)])
    sns.distplot(y1[~np.isnan(y1)], hist=True, kde=False, color="red",
                 kde_kws={'linewidth': 2},
                 label="|Corr| of k in TADs Mean: %.2f" % my1, ax=ax)
    my2 = np.mean(y2[~np.isnan(y2)])
    sns.distplot(y2[~np.isnan(y2)], hist=True, kde=False, color="blue",
                 kde_kws={'linewidth': 2},
                 label="|Corr| of k in others Mean: %.2f" % my2, ax=ax)
    ax.set_ylabel('#Regions')
    ax.set_xlabel('Genomic Distance(bp)')
    ax.legend()
    ax.set_xlim([0, 1])
    plt.savefig(fig_fp, dpi=300)


def correlation_between_k_in_TAD_bounaries_and_k_in_non_TAD(TAD_fp, non_TAD_fp):
    tad_df = pd.read_csv(TAD_fp, sep="\t", header=None).values
    non_tad_df = pd.read_csv(non_TAD_fp, sep="\t", header=None).values  # !!!!
    print("read successul")
    locis = {"chr%d" % chr_i: [] for chr_i in range(1, nchr + 1)}
    Ks = {"chr%d" % chr_i: [] for chr_i in range(1, nchr + 1)}

    n_samples_per_chromosome_for_control = 100  # number of regions as comparison

    for item in non_tad_df:
        chr_i, pos, val = item[0], item[1], item[4]
        locis[chr_i].append(int(pos))
        Ks[chr_i].append(val)

    n_tads = np.max(np.array(tad_df[:, -2]).astype(int))  # number of TADs
    print("%d TADs" % n_tads)
    d_corr = {}
    for i in range(1, n_tads):
        print("TAD %d/%d" % (i, n_tads))
        items = tad_df[tad_df[:, -2] == i]
        if items.shape[0] > 1:
            items_1 = items[items[:, -1] == 1]
            items_2 = items[items[:, -1] == 2]

            n_cpg1 = items_1.shape[0]
            n_cpg2 = items_2.shape[0]
            min_int_ncpg = min(n_cpg1, n_cpg2)
            max_int_ncpg = max(n_cpg1, n_cpg2)



            if min_int_ncpg > min_nCpG:
                # Calculating correlation of the CpGs in the TAD boundries
                corrs = []
                n_sites_avail_for_iter = max_int_ncpg - min_int_ncpg + 1
                distance = items_2[0][-3] - items_1[0][-3]

                ks1 = [items_1[ii][4] for ii in range(n_cpg1)]
                ks2 = [items_2[ii][4] for ii in range(n_cpg2)]
                fixed_ks = ks2 if n_cpg1 > n_cpg2 else ks1
                iter_ks = ks1 if n_cpg1 > n_cpg2 else ks2

                for ix in range(n_sites_avail_for_iter):
                    iks = iter_ks[ix: min_int_ncpg + ix]
                    ntc = abs(
                        np.corrcoef(fixed_ks, iks)[0][1])
                    corrs.append(ntc)
                corr_in_tads = np.mean(np.array(corrs))


                corr_samples = []

                sampling_rounds = 0
                for chr_i in range(1, 3):
                    ci = "chr%d" % chr_i
                    locis_chr = locis[ci]
                    Ks_chr = Ks[ci]
                    n_items = len(locis_chr)
                    sample_counter = 0
                    while sample_counter < n_samples_per_chromosome_for_control:
                        if sampling_rounds > 5:
                            break
                        sampling_rounds += 1
                        samples = random.sample(list(enumerate(locis_chr[0: int(n_items / 2)])),
                                                min(5000, int(n_items / 2)))
                        for ind_l1_r1, loci1_r1 in samples:
                            ind_l2_r1 = ind_l1_r1 + 1
                            loci2_r1 = locis_chr[ind_l2_r1]
                            cpg_num_r1 = 0
                            while loci2_r1 < loci1_r1 + window_sz * 2:
                                cpg_num_r1 += 1
                                ind_l2_r1 += 1
                                loci2_r1 = locis_chr[ind_l2_r1]

                            if cpg_num_r1 > min_nCpG:
                                loci1_r2 = loci1_r1 + distance
                                ind_l1_r2 = bisect_left(locis_chr, loci1_r2)
                                loci1_r2 = locis_chr[ind_l1_r2]
                                if abs(loci1_r2 - loci1_r1 - distance) < window_sz * 2:
                                    ind_l2_r2 = ind_l1_r2 + 1
                                    cpg_num_r2 = 0
                                    while locis_chr[ind_l2_r2] < loci1_r2 + window_sz * 2:
                                        cpg_num_r2 += 1
                                        ind_l2_r2 += 1
                                    if cpg_num_r2 > min_nCpG:
                                        ks_r1 = np.array(Ks_chr[ind_l1_r1: ind_l2_r1])
                                        ks_r2 = np.array(Ks_chr[ind_l1_r2: ind_l2_r2])
                                        l1 = ks_r1.shape[0]
                                        l2 = ks_r2.shape[0]
                                        min_int_ncpg = min(l1, l2)
                                        max_int_ncpg = max(l1, l2)

                                        corrs = []
                                        n_sites_avail_for_iter = max_int_ncpg - min_int_ncpg + 1

                                        fixed_ks = ks_r2 if l1 > l2 else ks_r1
                                        iter_ks = ks_r1 if l1 > l2 else ks_r2

                                        for ix in range(n_sites_avail_for_iter):
                                            iks = iter_ks[ix: min_int_ncpg + ix]
                                            ntc = abs(
                                                np.corrcoef(fixed_ks, iks)[0][1])
                                            corrs.append(ntc)
                                        non_tad_corr = np.mean(np.array(corrs))
                                        sample_counter += 1
                                        corr_samples.append(non_tad_corr)
                    if sampling_rounds > 5:
                        break
                if sampling_rounds < 5:
                    corr_for_non_tads = np.mean(np.array(corr_samples))
                    try:
                        if distance not in d_corr.keys():
                            d_corr[distance] = [corr_in_tads, corr_for_non_tads]
                        else:
                            tc1, nc1 = d_corr[distance]
                            d_corr[distance] = [(corr_in_tads + tc1) / 2., (corr_for_non_tads + nc1) / 2.]
                    except Exception as e:
                        pass
    dists = np.array(sorted(d_corr.keys()))
    tad_corrs = np.array([d_corr[d][0] for d in dists])
    non_tad_corrs = np.array([d_corr[d][1] for d in dists])
    return [dists, tad_corrs, non_tad_corrs]

def correlation_between_k_in_TAD_bounaries_and_k_in_non_TAD_large_slide_wd(TAD_fp, non_TAD_fp):
    tad_df = pd.read_csv(TAD_fp, sep="\t", header=None).values
    non_tad_df = pd.read_csv(non_TAD_fp, sep="\t", header=None).values  # !!!!
    print("read successul")
    locis = {"chr%d" % chr_i: [] for chr_i in range(1, nchr + 1)}
    Ks = {"chr%d" % chr_i: [] for chr_i in range(1, nchr + 1)}

    n_samples_per_chromosome_for_control = 100  # number of regions as comparison

    for item in non_tad_df:
        chr_i, pos, val = item[0], item[1], item[4]
        locis[chr_i].append(int(pos))
        Ks[chr_i].append(val)

    n_tads = np.max(np.array(tad_df[:, -2]).astype(int))  # number of TADs
    print("%d TADs" % n_tads)
    d_corr = {}
    for i in range(1, n_tads):
        print("TAD %d/%d" % (i, n_tads))
        items = tad_df[tad_df[:, -2] == i]
        if items.shape[0] > 1:
            items_1 = items[items[:, -1] == 1]
            items_2 = items[items[:, -1] == 2]

            n_cpg1 = items_1.shape[0]
            n_cpg2 = items_2.shape[0]
            min_int_ncpg = min(n_cpg1, n_cpg2)
            max_int_ncpg = max(n_cpg1, n_cpg2)

            if min_int_ncpg > min_nCpG:
                # Calculating correlation of the CpGs in the TAD boundries
                corrs = []
                n_sites_avail_for_iter = max_int_ncpg - min_int_ncpg + 1
                distance = items_2[0][-3] - items_1[0][-3]

                ks1 = [items_1[ii][4] for ii in range(n_cpg1)]
                ks2 = [items_2[ii][4] for ii in range(n_cpg2)]
                fixed_ks = ks2 if n_cpg1 > n_cpg2 else ks1
                iter_ks = ks1 if n_cpg1 > n_cpg2 else ks2

                for ix in range(n_sites_avail_for_iter):
                    iks = iter_ks[ix: min_int_ncpg + ix]
                    ntc = abs(
                        np.corrcoef(fixed_ks, iks)[0][1])
                    corrs.append(ntc)
                corr_in_tads = np.mean(np.array(corrs))


                corr_samples = []

                sampling_rounds = 0
                for chr_i in range(1, 3):
                    ci = "chr%d" % chr_i
                    locis_chr = locis[ci]
                    Ks_chr = Ks[ci]
                    n_items = len(locis_chr)
                    sample_counter = 0
                    while sample_counter < n_samples_per_chromosome_for_control:
                        if sampling_rounds > 5:
                            break
                        sampling_rounds += 1
                        samples = random.sample(list(enumerate(locis_chr[0: int(n_items / 2)])),
                                                min(5000, int(n_items / 2)))
                        for ind_l1_r1, loci1_r1 in samples:
                            ind_l2_r1 = ind_l1_r1 + 1
                            loci2_r1 = locis_chr[ind_l2_r1]
                            cpg_num_r1 = 0
                            while loci2_r1 < loci1_r1 + window_sz * 2:
                                cpg_num_r1 += 1
                                ind_l2_r1 += 1
                                loci2_r1 = locis_chr[ind_l2_r1]

                            if cpg_num_r1 > min_nCpG:
                                loci1_r2 = loci1_r1 + distance
                                ind_l1_r2 = bisect_left(locis_chr, loci1_r2)
                                loci1_r2 = locis_chr[ind_l1_r2]
                                if abs(loci1_r2 - loci1_r1 - distance) < window_sz * 2:
                                    ind_l2_r2 = ind_l1_r2 + 1
                                    cpg_num_r2 = 0
                                    while locis_chr[ind_l2_r2] < loci1_r2 + window_sz * 2:
                                        cpg_num_r2 += 1
                                        ind_l2_r2 += 1
                                    if cpg_num_r2 > min_nCpG:
                                        ks_r1 = np.array(Ks_chr[ind_l1_r1: ind_l2_r1])
                                        ks_r2 = np.array(Ks_chr[ind_l1_r2: ind_l2_r2])
                                        l1 = ks_r1.shape[0]
                                        l2 = ks_r2.shape[0]
                                        min_int_ncpg = min(l1, l2)
                                        max_int_ncpg = max(l1, l2)

                                        corrs = []
                                        n_sites_avail_for_iter = max_int_ncpg - min_int_ncpg + 1

                                        fixed_ks = ks_r2 if l1 > l2 else ks_r1
                                        iter_ks = ks_r1 if l1 > l2 else ks_r2

                                        for ix in range(n_sites_avail_for_iter):
                                            iks = iter_ks[ix: min_int_ncpg + ix]
                                            ntc = abs(
                                                np.corrcoef(fixed_ks, iks)[0][1])
                                            corrs.append(ntc)
                                        non_tad_corr = np.mean(np.array(corrs))
                                        sample_counter += 1
                                        corr_samples.append(non_tad_corr)
                    if sampling_rounds > 5:
                        break
                if sampling_rounds < 5:
                    corr_for_non_tads = np.mean(np.array(corr_samples))
                    try:
                        if distance not in d_corr.keys():
                            d_corr[distance] = [corr_in_tads, corr_for_non_tads]
                        else:
                            tc1, nc1 = d_corr[distance]
                            d_corr[distance] = [(corr_in_tads + tc1) / 2., (corr_for_non_tads + nc1) / 2.]
                    except Exception as e:
                        pass
    dists = np.array(sorted(d_corr.keys()))
    tad_corrs = np.array([d_corr[d][0] for d in dists])
    non_tad_corrs = np.array([d_corr[d][1] for d in dists])
    return [dists, tad_corrs, non_tad_corrs]


if __name__ == "__main__":
    TAD_fp = os.path.join(data_dir, "TADs", "K_inter_TADs_%dbp.bed" %  window_sz)
    non_TAD_fp = os.path.join(data_dir, "TADs", "K_non_TADs.bed")
    histogram_of_TAD_length(TAD_fp, non_TAD_fp)


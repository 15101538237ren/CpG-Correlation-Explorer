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
D_MAX = 1000
BASE_DIR = ".."
ChrmoHMM_LABELS =["Active Promoter", "Weak Promoter", "Poised Promoter", "Strong Enhancer", "Strong Enhancer", "Weak Enhancer", "Weak Enhancer", "Insulator", "Txn Transition", "Txn Elongation", "Weak Txn", "Repressed", "Heterochrom_lo", "Repetitive_CNV", "Repetitive_CNV"]
FILE_ORDERED_NAMES = {
            "ChromHMM": ["1_Active_Promoter", "2_Weak_Promoter", "3_Poised_Promoter", "4_Strong_Enhancer", "5_Strong_Enhancer",
                         "6_Weak_Enhancer", "7_Weak_Enhancer", "8_Insulator", "9_Txn_Transition", "10_Txn_Elongation",
                         "11_Weak_Txn", "12_Repressed", "13_Heterochrom_lo", "14_Repetitive_CNV", "15_Repetitive_CNV"],
            "Genomic_Regions": ["Promoter", "Enhancer", "CGI", "Exons", "Introns", "Intergenic", "LINE", "SINE", "LTR", "5UTR", "3UTR",],
            "Histone_Modification": ["H3k4me1Not3_no_k27me3_no_k27ac", "H3k4me1Not3_k27me3", "H3k4me1Not3_k27ac",
                                                    "H3k4me3Not1_no_k27me3_no_k27ac", "H3k4me3Not1_k27me3", "H3k4me3Not1_k27ac",
                                                    "H3K9Me3", "H3K36me3"],
            "TFBS" : ["DNase", "EZH2", "H2AZ", "tfbs_cluster_v3"]
            }
FILE_LABELS = {
            "ChromHMM": ChrmoHMM_LABELS,
            "Genomic_Regions": ["Promoter", "Enhancer", "CGI", "Exons", "Introns", "Intergenic", "LINE", "SINE", "LTR", "5UTR", "3UTR"],
            "Histone_Modification": ["H3k4me1", "H3k4me1 + H3k27me3", "H3k4me1 + H3k27ac",
                                        "H3k4me3", "H3k4me3 + H3k27me3", "H3k4me3 + H3k27ac", "H3K9Me3", "H3K36me3"],
            "TFBS": ["DNase", "EZH2", "H2AZ", "TFBS"]
                }

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def plot_whole_landscape(fig_format="png"):
    regions = ["Genomic_Regions", "Histone_Modification", "ChromHMM", "TFBS"]  #
    RD_DIRNAME = "Whole_Landscape"
    cm = plt.get_cmap('gist_rainbow')
    FIG_DIR = os.path.join(BASE_DIR, "figures")
    NBIN = 30
    vmins = [0, 120]
    vmaxs = [0.15, 220]
    N_COL = 18
    COL_Labels = ["Corr", "Corr with DNase", "Corr with Nuc Occ", "Hist of k", "Hist of methylation", "Hist2D of k/methy",
                  "DNase Peak", "Nucleosome Occupancy", "H3k4me1", "H3k4me3", "H3k9me3", "H3k9ac", "H3k27ac", "H3k27me3",
                  "H3k36me3", "H4k20me", "CTCF", "p300"]
    COL_INDEXS = [0, 0, 0, 4, 3, 0, -2, -1, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    MP = 50 #Max Peak Value
    xlims = [D_MAX, D_MAX, D_MAX, 10, 1., 1, 0.15, 220, MP, MP, MP, MP, MP, MP, MP, MP, MP, MP]
    ylims = [1, 1, 1, 0.6, 16, 1, 70, 0.05, 0.25, 0.2, 0.14, 0.06, .1, .12, .25, .20, 0.05, .2]
    for REGION in regions:
        fig_dir = os.path.join(FIG_DIR, RD_DIRNAME)
        mkdir(fig_dir)
        out_rd_corr_dir = os.path.join("../data/K_region_intersect", REGION, "K_Rd")
        file_paths = [os.path.join(out_rd_corr_dir, "%s.bed" % file_name) for file_name in FILE_ORDERED_NAMES[REGION]]
        file_labels = FILE_LABELS[REGION]
        bed_fps = [os.path.join("../data/K_region_intersect", REGION, "%s.bed" % file_name) for file_name in FILE_ORDERED_NAMES[REGION]]
        N_ROW = len(file_labels)

        fig_fp = os.path.join(fig_dir, "%s.%s" % (REGION,  fig_format))
        fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
        for row in range(N_ROW):
            RD_df = pd.read_csv(file_paths[row], sep="\t", header=None).values
            bed_df = pd.read_csv(bed_fps[row], sep="\t", header=None).values
            for col in range(N_COL):
                print("%s, %s" % (file_labels[row], COL_Labels[col]))
                ax = axs[row][col]
                xlim = xlims[col]
                ylim = ylims[col]
                if col <= 2:
                    x = RD_df[:, 0]
                    y = RD_df[:, 1]
                    try:
                        y2 = localreg(x, y, degree=2, kernel=tricube, width=100)
                        if col != 0:
                            z = RD_df[:, 1 + col]
                            sc = ax.scatter(x, y, s=8, c=z, label=file_labels[row], cmap=cm, vmin=vmins[col - 1],
                                            vmax=vmaxs[col - 1])
                            fig.colorbar(sc, ax=ax)
                        else:
                            ax.scatter(x, y, s=2, color=cm(1. * row / N_ROW), label=file_labels[row])
                        ax.plot(x, y2, "k-", linewidth=2)
                    except np.linalg.LinAlgError as e:
                        sns.regplot(x=x, y=y, ax=ax, scatter_kws={'s': 8, 'color': cm(1. * row / N_ROW)})  # , line_kws ={'color':'black', "lw": 2}
                    ax.set_xticks(range(0, D_MAX + 1, 200))
                    ax.set_xlim(0, xlim)
                    ax.set_ylim(0, ylim)
                    ax.set_title(COL_Labels[col], fontsize=18)
                    if col == 0:
                        ax.set_ylabel(file_labels[row], fontsize=20)
                    if row == N_ROW - 1:
                        ax.set_xlabel("Genomic Distance(bp)")
                else:
                    col_ind_in_bed = COL_INDEXS[col]
                    if col_ind_in_bed != 0:
                        vals = bed_df[:, col_ind_in_bed]
                        vals = vals[vals != "."].astype(float)
                        vals[vals > xlim] = xlim
                        _ = ax.hist(vals, bins=NBIN, density=True, color=cm(1. * col / N_COL), edgecolor='black',
                                    alpha=0.5, linewidth=0.5)
                        ax.set_xlim([0, xlim])
                        ax.set_ylim(0, ylim)
                        ax.set_title(COL_Labels[col], fontsize=18)
                    else:
                        ks = bed_df[:, 4].astype(float)
                        methys = bed_df[:, 3].astype(float)
                        ks[ks <= 0] = 0.001
                        ks = np.log10(ks)
                        h = ax.hist2d(methys, ks, bins=(NBIN, NBIN), density=True, vmin=0, vmax=5, cmap="viridis")

                        fig.colorbar(h[3], ax=ax)
                        if col == 0:
                            ax.set_ylabel('log10(K)')
                        if row == N_ROW - 1:
                            ax.set_xlabel('Methylation level')
                        ax.set_ylim(-1, 1)
                        ax.set_xlim(0, 1)
                        ax.set_xticks([0.2 * i for i in range(6)])
        plt.savefig(fig_fp, dpi=300, bbox_inches='tight', pad_inches=0.1)
if __name__ == "__main__":
    plot_whole_landscape()
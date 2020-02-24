# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from localreg import *
import seaborn as sns

sns.set(color_codes=True)
data_dir = "../data"
figure_dir = "../figures/Q1_CGI_threshold"
EACH_SUB_FIG_SIZE = 5
D_MAX = 5000
RATIOS = ["06", "07", "08", "09", "10", "11", "12", "13"]
RATIO_LABELS = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def partition_bed_by_region_labels_and_generate_config_file(config_fp, K_RD= 1):
    correlation_type = '' if K_RD else '-methy'

    out_dir = os.path.join(data_dir, "CGI_identified_with_different_thereshold")
    mkdir(out_dir)
    with open(config_fp, "w") as config:
        for ratio in RATIOS:
            input_fp = os.path.join(out_dir, "CGI_%s_K_intersected" % ratio + correlation_type +".bed")
            Rd_fp = os.path.join(out_dir, "CGI_%s_K_intersected" % ratio + correlation_type +"-only-within-Rd.bed" )
            config.write("%s\t%s\t%d\n" %(input_fp, Rd_fp, K_RD))
        print("%s saved" % config_fp)

def plot_local_regression_and_RD(fig_format="png"):
    K_RD = [1]# , 0
    N_COL = 4
    N_ROW = 2
    cm = plt.get_cmap('gist_rainbow')
    mkdir(figure_dir)
    for km in K_RD:
        correlation_type = '' if km else '-methy'
        fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
        fig_fp = os.path.join(figure_dir, "CGI_RATIO_FROM-06-13-%d.%s" % (D_MAX, fig_format))
        for rid, ratio in enumerate(RATIOS):
            correlation_fp = os.path.join(data_dir, "CGI_identified_with_different_thereshold", "CGI_%s_K_intersected" % ratio + correlation_type +"-only-within-Rd.bed")
            file_label = "Obs/Expect " + str(RATIO_LABELS[rid])
            row = rid // N_COL
            col = rid % N_COL
            ax = axs[row][col]
            RD_df = pd.read_csv(correlation_fp, sep="\t", header=None).values
            x = RD_df[:, 0]
            y = RD_df[:, 1]
            color = cm(1. * rid / len(RATIOS))
            try:
                y2 = localreg(x, y, degree=2, kernel=tricube, width=100)
                ax.scatter(x, y, s=8, color=color, label=file_label)
                ax.plot(x, y2, "k-", linewidth=2)
            except np.linalg.LinAlgError as e:
                sns.regplot(x=x, y=y, ax=ax, scatter_kws={'s':8, 'color': color})#, line_kws ={'color':'black', "lw": 2}
            ax.set_xticks(range(0, D_MAX + 1, 1000))
            ax.set_xlim(0, D_MAX)
            ax.set_ylim(0, 1.0)
            ax.set_title(file_label, fontsize=14)
        plt.savefig(fig_fp, dpi=300, bbox_inches='tight', pad_inches=0.1)

def plot_local_regression_and_RD_within_vs_across(fig_format="png"):
    K_RD = [1]# , 0
    N_COL = 4
    N_ROW = 2
    cm = plt.get_cmap('gist_rainbow')
    mkdir(figure_dir)
    for km in K_RD:
        correlation_type = '' if km else '-methy'
        fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
        fig_fp = os.path.join(figure_dir, "CGI_RATIO_Within-Across-5000bp.%s" % fig_format)
        for rid, ratio in enumerate(RATIOS):
            within_correlation_fp = os.path.join(data_dir, "CGI_identified_with_different_thereshold", "CGI_%s_K_intersected" % ratio + correlation_type +"-only-within-Rd.bed")
            across_correlation_fp = os.path.join(data_dir, "CGI_identified_with_different_thereshold", "CGI_%s_K_intersected" % ratio + correlation_type +"-both-within-across-Rd.bed")
            file_label = "Obs/Expect " + str(RATIO_LABELS[rid])
            row = rid // N_COL
            col = rid % N_COL
            ax = axs[row][col]
            within_RD = pd.read_csv(within_correlation_fp, sep="\t", header=None).values
            across_RD = pd.read_csv(across_correlation_fp, sep="\t", header=None).values
            x1 = within_RD[:, 0]
            y1 = within_RD[:, 1]

            x2 = across_RD[:, 0]
            y2 = across_RD[:, 1]

            color = cm(1. * rid / len(RATIOS))
            try:
                yy1 = localreg(x1, y1, degree=2, kernel=tricube, width=100)
                yy2 = localreg(x2, y2, degree=2, kernel=tricube, width=100)
                ax.scatter(x1, y1, s=3, color=color, label="within")
                ax.scatter(x2, y2, s=1.5, color="black", label="within")
                ax.plot(x2, yy2, "w-", linewidth=1.5)
                ax.plot(x1, yy1, "k-", linewidth=1.5)
            except np.linalg.LinAlgError as e:
                sns.regplot(x=x1, y=y1, ax=ax, scatter_kws={'s':8, 'color': color})#, line_kws ={'color':'black', "lw": 2}
                sns.regplot(x=x2, y=y2, ax=ax, scatter_kws={'s':8, 'color': "black"})#, line_kws ={'color':'black', "lw": 2}
            ax.set_xticks(range(0, D_MAX + 1, 1000))
            ax.set_xlim(0, D_MAX)
            ax.set_ylim(0, 1.0)
            ax.set_title(file_label, fontsize=14)
        plt.savefig(fig_fp, dpi=300, bbox_inches='tight', pad_inches=0.1)

def plot_hist_of_k_methy_or_len_CGI(property_name, plotting_together= True, fig_format="png"):

    NBIN = 30
    cm = plt.get_cmap('gist_rainbow')
    mkdir(figure_dir)
    fig_fp = os.path.join(figure_dir, "%s_FROM-06-13.%s" % (property_name, fig_format))
    if plotting_together:
        N_COL = 1
        N_ROW = 1
    else:
        N_COL = 4
        N_ROW = 2
    fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
    if property_name == "k":
        col_idx = 4
        xlim = 10.
        ylim = 1.3
    elif property_name == "methy":
        col_idx = 3
        xlim = 1.
        ylim = 20
    else:
        col_idx = 3
        xlim = 5000.
        ylim = 0.0020
    for rid, ratio in enumerate(RATIOS):
        if property_name != "len-CGI":
            bed_fp = os.path.join(data_dir, "CGI_identified_with_different_thereshold", "CGI_%s_K_intersected" % ratio + ".bed")
        else:
            bed_fp = os.path.join(data_dir, "CGI_identified_with_different_thereshold",
                                  "CGI_%s_ratio" % ratio + ".bed")
        file_label = "Obs/Exp " + str(RATIO_LABELS[rid])
        row = rid // N_COL
        col = rid % N_COL

        df = pd.read_csv(bed_fp, sep="\t", header=None).values
        vals = df[:, col_idx]

        color = cm(1. * rid / len(RATIOS))

        if plotting_together:
            ax = axs
            sns.distplot(vals, ax=ax, hist=False, kde=True,
                         kde_kws={'linewidth': 2}, color= color,
                         label=file_label)
        else:
            ax = axs[row][col]
            ax.hist(vals, bins=NBIN, range=(0, xlim), density=True, color=color, edgecolor='black', alpha=0.5, linewidth=0.5)

        ax.set_xticks([round(0. + i*(xlim/5), 2) for i in range(5)])
        if row == N_ROW - 1:
            ax.set_xlabel(property_name)
        ax.set_xlim(0, xlim)
        ax.set_ylim(0, ylim)
        ax.set_title(file_label, fontsize=14)
    plt.savefig(fig_fp, dpi=300, bbox_inches='tight', pad_inches=0.1)

def plot_heatmap_of_k_methy_for_diff_CGI_threshold(fig_format="png"):

    NBIN = 30
    cm = plt.get_cmap('gist_rainbow')
    mkdir(figure_dir)
    fig_fp = os.path.join(figure_dir, "heatmap.%s" % (fig_format))
    N_COL = 4
    N_ROW = 2
    fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))

    for rid, ratio in enumerate(RATIOS):
        bed_fp = os.path.join(data_dir, "CGI_identified_with_different_thereshold",
                                  "CGI_%s_K_intersected" % ratio + ".bed")
        file_label = "Obs/Exp " + str(RATIO_LABELS[rid])
        row = rid // N_COL
        col = rid % N_COL
        ax = axs[row][col]

        df = pd.read_csv(bed_fp, sep="\t", header=None).values
        methys = df[:, 3].astype(float)
        ks = df[:, 4].astype(float)
        ks[ks <= 0] = 0.001
        ks = np.log10(ks)

        corr_coef = np.corrcoef(ks, methys)[0][1]
        h = ax.hist2d(methys, ks, bins=(NBIN, NBIN), density=True, vmin=0, vmax=5, cmap="viridis")

        fig.colorbar(h[3], ax=ax)
        if col == 0:
            ax.set_ylabel('log10(K)')
        if row == N_ROW - 1:
            ax.set_xlabel('Methylation level')
        ax.set_ylim(-1, 1)
        ax.set_xlim(0, 1)

        ax.set_xticks([0.2 *i for i in range(6)])
        ax.set_title(file_label + ", Corr:%.2f" % corr_coef, fontsize=14)
    plt.savefig(fig_fp, dpi=300, bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":
    plot_local_regression_and_RD()
    # plot_heatmap_of_k_methy_for_diff_CGI_threshold()
    # config_fp = os.path.join(data_dir, "cgi_threshold-only-within.config")
    # partition_bed_by_region_labels_and_generate_config_file(config_fp)
    # for property_name in ["len-CGI", "k", "methy"]:#
    #     plot_hist_of_k_methy_or_len_CGI(property_name, plotting_together=True)

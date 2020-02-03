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
D_MAX =1000
RATIO_LABELS = ["0.%d" % i for i in range(3, 10)]

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def partition_bed_by_region_labels_and_generate_config_file(config_fp, K_RD= 1):
    correlation_type = '' if K_RD else '-methy'

    out_dir = os.path.join(data_dir, "CGI_identified_with_different_thereshold")
    mkdir(out_dir)
    with open(config_fp, "w") as config:
        for ratio in range(3, 10):
            input_fp = os.path.join(out_dir, "CGI_0%d_K_intersected" % ratio + correlation_type +".bed")
            Rd_fp = os.path.join(out_dir, "CGI_0%d_K_intersected" % ratio + correlation_type +"-Rd.bed" )
            config.write("%s\t%s\t%d\n" %(input_fp, Rd_fp, K_RD))
        print("%s saved" % config_fp)

def plot_local_regression_and_RD(fig_format="png"):
    K_RD = [1]# , 0
    N_COL = 4
    N_ROW = 2
    RATIOS = range(3, 10)
    cm = plt.get_cmap('gist_rainbow')
    mkdir(figure_dir)
    for km in K_RD:
        correlation_type = '' if km else '-methy'
        fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
        fig_fp = os.path.join(figure_dir, "CGI_RATIO_FROM-03-09.%s" % fig_format)
        for rid, ratio in enumerate(RATIOS):
            correlation_fp = os.path.join(data_dir, "CGI_identified_with_different_thereshold", "CGI_0%d_K_intersected" % ratio + correlation_type +"-Rd.bed")
            file_label = "Obs/Expect " + RATIO_LABELS[rid]
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
            ax.set_xticks(range(0, D_MAX + 1, 200))
            ax.set_xlim(0, D_MAX)
            ax.set_ylim(0, 1.0)
            ax.set_title(file_label, fontsize=14)
        plt.savefig(fig_fp, dpi=300, bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":
    config_fp = os.path.join(data_dir, "cgi_threshold.config")
    #partition_bed_by_region_labels_and_generate_config_file(config_fp)
    plot_local_regression_and_RD()

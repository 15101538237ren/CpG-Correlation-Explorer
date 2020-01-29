# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from localreg import *
import seaborn as sns

sns.set(color_codes=True)
data_dir = "../data"
figure_dir = "../figures/Q1_CGI_density"
EACH_SUB_FIG_SIZE = 5
D_MAX =1000
REGION_LABELS = ["Promoter", "Enhancer", "f-UTR"]
CGI_LABELS = ["Non-CGI", "CGI"]

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def partition_bed_by_region_labels_and_generate_config_file(config_fp, K_RD= 1):
    correlation_type = '' if K_RD else '-methy'
    input_bed_fp = os.path.join(data_dir, "K_CGI_Prom_Enh_5UTR.bed")
    df = pd.read_csv(input_bed_fp, sep='\t', header=None).values
    with open(config_fp, "w") as config:
        for rid, region_label in enumerate(REGION_LABELS):
            out_dir = os.path.join(data_dir, "Q1_CGI_density", region_label)
            mkdir(out_dir)

            REGIONS = ["Non-" + region_label, region_label]

            for iid, region in enumerate(REGIONS):
                for jid, cgi in enumerate(CGI_LABELS):
                    cgi_col_idx = 6
                    region_col_idx = cgi_col_idx + rid + 1
                    out_fp = os.path.join(out_dir, region + "-" + cgi + correlation_type +".bed")
                    Rd_fp = os.path.join(out_dir, region + "-" + cgi + correlation_type +"-Rd.bed")
                    target_rows = np.logical_and(df[:, cgi_col_idx].astype(int) == jid, df[:, region_col_idx].astype(int) == iid)

                    target_df = df[target_rows]
                    np.savetxt(out_fp, target_df[:, 0: cgi_col_idx], fmt='%s\t%d\t%d\t%f\t%f\t%f', delimiter='\n')
                    print("%s saved" % out_fp)
                    config.write("%s\t%s\t%d\n" %(out_fp, Rd_fp, K_RD))

def plot_local_regression_and_RD(fig_format="png"):
    K_RD = [1]# , 0
    N_COL = N_ROW =  2
    cm = plt.get_cmap('gist_rainbow')
    mkdir(figure_dir)
    for km in K_RD:
        correlation_type = '' if km else '-methy'
        for rid, region_label in enumerate(REGION_LABELS):
            input_dir = os.path.join(data_dir, "Q1_CGI_density", region_label)

            REGIONS = ["Non-" + region_label, region_label]

            fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
            fig_fp = os.path.join(figure_dir, "%s.%s" % (region_label, fig_format))
            for iid, region in enumerate(REGIONS):
                for jid, cgi in enumerate(CGI_LABELS):
                    correlation_fp = os.path.join(input_dir, region + "-" + cgi + correlation_type +"-Rd.bed")
                    file_label = region + " " + cgi
                    idx = iid * len(REGIONS) + jid
                    row = iid
                    col = jid
                    ax = axs[row][col]
                    RD_df = pd.read_csv(correlation_fp, sep="\t", header=None).values
                    x = RD_df[:, 0]
                    y = RD_df[:, 1]
                    color = cm(1. * idx / (N_COL*N_ROW))
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
    config_fp = os.path.join(data_dir, "cgi_non_cgi.config")
    partition_bed_by_region_labels_and_generate_config_file(config_fp)
    #plot_local_regression_and_RD()

# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from localreg import *
import seaborn as sns

sns.set(color_codes=True)
data_dir = "../data"
figure_dir = "../figures/nucleosome"
EACH_SUB_FIG_SIZE = 5
D_MAX = 3000


def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def plot_nucleosome_occupancy_sample(fig_format ="png"):
    mkdir(figure_dir)
    input_fp = os.path.join(data_dir, "nucleosome_positioning", "GSM1194220_H1.bed")
    df = pd.read_csv(input_fp, sep="\t", header=None).values
    NROW = 4
    NCOL = 4
    fig, axs = plt.subplots(NROW, NROW, figsize=(NROW * EACH_SUB_FIG_SIZE, NROW * EACH_SUB_FIG_SIZE))
    fig_fp = os.path.join(figure_dir, "nucleosome_occupancy_of_part_of_chr1.%s" % fig_format)
    for rid in range(NROW * NCOL):
        row = rid // NCOL
        col = rid % NCOL
        ax = axs[row][col]
        start = 20000 + rid * 2000
        end = 20000 + (rid + 1) * 2000
        loci = df[start:end, 1].astype(np.float)
        oc = df[start:end, 2].astype(np.float)
        try:
            y2 = localreg(loci, oc, degree=2, kernel=tricube, width=10)
            ax.scatter(loci, oc, s=8, color="blue")
            ax.plot(loci, y2, "k-", linewidth=2)
        except np.linalg.LinAlgError as e:
            sns.regplot(x=loci, y=oc, ax=ax, scatter_kws={'s': 8, 'color': "blue"})
        ax.set_xlim(start, end)
        ax.set_ylim(0, 600)
        ax.set_title("nucleos occupancy of chr1: %d:%d" %(start, end) , fontsize=14)


    plt.savefig(fig_fp, dpi=300, bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":
    plot_nucleosome_occupancy_sample()
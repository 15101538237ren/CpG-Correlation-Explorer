# -*- coding: utf-8 -*-
import os, pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from localreg import *

data_dir = "data"
chromHMM_labels = range(1, 16)#
DMAX = 1000
distance_range = range(2, DMAX + 1)
chromosomes = ["chr%d" % i for i in range(1, 2)]
FEATURE_NAME = "feature"
CORR_NAME = "corr"
figure_dir = "figures/kdiff_vs_feature"
EACH_SUB_FIG_SIZE = 5
ChrmoHMM_LABELS =["Active Promoter", "Weak Promoter", "Poised Promoter", "Strong Enhancer", "Strong Enhancer", "Weak Enhancer", "Weak Enhancer", "Insulator", "Txn Transition", "Txn Elongation", "Weak Txn", "Repressed", "Heterochrom_lo", "Repetitive_CNV", "Repetitive_CNV"]

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def plot_corr_vs_distance(pkl_fp):
    mkdir(figure_dir)
    cm = plt.get_cmap('gist_rainbow')
    N_COL = 5
    N_ROW = 3
    with open(pkl_fp, 'rb') as pkl:
        data_list = pickle.load(pkl)
        fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))

        for j, chrHMMlbl in enumerate(chromHMM_labels):
            row = j // N_COL
            col = j % N_COL
            if N_ROW == 1:
                ax = axs[col]
            else:
                ax = axs[row][col]
            corr = []
            features = []
            for dataset in data_list:
                dataset = dataset[0]
                for d in distance_range:
                    corr.append(dataset[chrHMMlbl][d][CORR_NAME])
                    features.append(d)
            ax.scatter(features, corr, color=cm(1. * chrHMMlbl / 15.), s=3)
            ax.set_title(ChrmoHMM_LABELS[j])
            ax.set_xlim([0, DMAX])
            ax.set_ylim([-1, 1.0])
            if row == N_ROW - 1:
                ax.set_xlabel("1D-distance", fontsize= 12)
            if col == 0:
                ax.set_ylabel("Correlation", fontsize= 12)
        plt.savefig(os.path.join(figure_dir, "distance.png"), dpi=300)

def plot_corr_vs_local_CpG_density(pkl_fp, split= False):
    mkdir(figure_dir)
    cm = plt.get_cmap('gist_rainbow')
    N_COL = 5
    N_ROW = 3
    d = 10
    with open(pkl_fp, 'rb') as pkl:
        data_list = pickle.load(pkl)
        if split:
            fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
            for j, chrHMMlbl in enumerate(chromHMM_labels):
                row = j // N_COL
                col = j % N_COL
                if N_ROW == 1:
                    ax = axs[col]
                else:
                    ax = axs[row][col]
                corr = []
                features = []
                for dataset in data_list:
                    dataset = dataset[0]
                    feature_data = dataset[chrHMMlbl][d][FEATURE_NAME]

                    if feature_data:
                        for feature in feature_data:
                            features.append((feature[0][-2] + feature[1][-2]) * 0.5)
                            corr.append(dataset[chrHMMlbl][d][CORR_NAME])
                ax.scatter(features, corr, color=cm(1. * chrHMMlbl / 15.), s=3)
                ax.set_title(ChrmoHMM_LABELS[j])
                ax.set_xlim([0, 1.0])
                ax.set_ylim([-1, 1.0])
                if row == N_ROW - 1:
                    ax.set_xlabel("Local CpG Density", fontsize= 12)
                if col == 0:
                    ax.set_ylabel("Correlation", fontsize= 12)
            plt.savefig(os.path.join(figure_dir, "local_CpG_density_for_d-%d_split.png" % d), dpi=300)
        else:
            fig, ax = plt.subplots(1, 1, figsize=(1 * EACH_SUB_FIG_SIZE, 1* EACH_SUB_FIG_SIZE))

            for i, chrHMMlbl in enumerate(chromHMM_labels):
                corr = []
                features = []
                for dataset in data_list:
                    dataset = dataset[0]
                    feature_data = dataset[chrHMMlbl][d][FEATURE_NAME]
                    if feature_data:
                        for feature in feature_data:
                            features.append((feature[0][-2] + feature[1][-2]) * 0.5)
                            corr.append(dataset[chrHMMlbl][d][CORR_NAME])
                ax.scatter(features, corr, color=cm(1 * i / 15.), s=3)
            ax.set_xlim([0, 1.0])
            ax.set_ylim([-1, 1.0])
            ax.set_xlabel("Local CpG Density", fontsize= 12)
            ax.set_ylabel("Correlation", fontsize= 12)
            plt.savefig(os.path.join(figure_dir, "local_CpG_density_for_d-%d.png" % d), dpi=300)

def plot_kdiff_vs_distance(pkl_fp):
    mkdir(figure_dir)
    cm = plt.get_cmap('gist_rainbow')
    N_COL = 5
    N_ROW = 3
    with open(pkl_fp, 'rb') as pkl:
        CpG_pair_withd_chrHMM = pickle.load(pkl)
        CpG_pair_withd_chrHMM = CpG_pair_withd_chrHMM[0]
        fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))

        for j, chrHMMlbl in enumerate(chromHMM_labels):
            row = j // N_COL
            col = j % N_COL
            if N_ROW == 1:
                ax = axs[col]
            else:
                ax = axs[row][col]
            kdiffs = []
            features = []
            kdiff_std = []
            for d in distance_range:
                feature_items = CpG_pair_withd_chrHMM[chrHMMlbl][d][FEATURE_NAME]
                kdiffs.append(np.mean(np.array([item[2] for item in feature_items])))
                kdiff_std.append(np.std(np.array([item[2] for item in feature_items])))
                features.append(d)

            ax.scatter(features, kdiffs, color=cm(1. * chrHMMlbl / 15.), s=3)
            #y2 = localreg(np.array(features), np.array(kdiff_std), degree=2, kernel=tricube, width=100)
            ax.scatter(features, kdiff_std, color='black', s=3)
            #ax.plot(np.array(features), y2, color='white')
            ax.set_title(ChrmoHMM_LABELS[j])
            ax.set_xlim([0, DMAX])
            ax.set_ylim([0, 8.0])
            if row == N_ROW - 1:
                ax.set_xlabel("1D-distance", fontsize= 12)
            if col == 0:
                ax.set_ylabel("abs(K1 - K2)", fontsize= 12)
        plt.savefig(os.path.join(figure_dir, "k_diff_distance.png"), dpi=300)

def plot_kdiff_vs_local_CpG_density(pkl_fp, split= False):
    mkdir(figure_dir)
    cm = plt.get_cmap('gist_rainbow')
    N_COL = 5
    N_ROW = 3
    d = 100
    with open(pkl_fp, 'rb') as pkl:
        CpG_pair_withd_chrHMM = pickle.load(pkl)
        CpG_pair_withd_chrHMM = CpG_pair_withd_chrHMM[0]
        if split:
            fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
            for j, chrHMMlbl in enumerate(chromHMM_labels):
                row = j // N_COL
                col = j % N_COL
                if N_ROW == 1:
                    ax = axs[col]
                else:
                    ax = axs[row][col]
                kdiffs = []
                features = []
                feature_data = CpG_pair_withd_chrHMM[chrHMMlbl][d][FEATURE_NAME]

                if feature_data:
                    for feature in feature_data:
                        features.append((feature[0][-2] + feature[1][-2]) * 0.5)
                        kdiffs.append(feature[2])
                ax.scatter(features, kdiffs, color=cm(1. * chrHMMlbl / 15.), s=3)
                ax.set_title(ChrmoHMM_LABELS[j])
                ax.set_xlim([0, 1.0])
                ax.set_ylim([0, 8.0])
                if row == N_ROW - 1:
                    ax.set_xlabel("Local CpG Density", fontsize= 12)
                if col == 0:
                    ax.set_ylabel("abs(K1 - K2)", fontsize= 12)
            plt.savefig(os.path.join(figure_dir, "local_CpG_density_for_d-%d_split.png" % d), dpi=300)
        else:
            fig, ax = plt.subplots(1, 1, figsize=(1 * EACH_SUB_FIG_SIZE, 1* EACH_SUB_FIG_SIZE))

            for i, chrHMMlbl in enumerate(chromHMM_labels):
                kdiffs = []
                features = []
                feature_data = CpG_pair_withd_chrHMM[chrHMMlbl][d][FEATURE_NAME]
                if feature_data:
                    for feature in feature_data:
                        features.append((feature[0][-2] + feature[1][-2]) * 0.5)
                        kdiffs.append(feature[2])
                ax.scatter(features, kdiffs, color=cm(1 * i / 15.), s=3)
            ax.set_xlim([0, 1.0])
            ax.set_ylim([0, 8.0])
            ax.set_xlabel("Local CpG Density", fontsize= 12)
            ax.set_ylabel("abs(K1 - K2)", fontsize= 12)
            plt.savefig(os.path.join(figure_dir, "local_CpG_density_with_kdiff_for_d-%d.png" % d), dpi=300)

if __name__ == "__main__":
    pkl_fp = os.path.join(data_dir, "CpG_pairs_with_d_ChromHMM.pkl")
    #plot_kdiff_vs_distance(pkl_fp)
    plot_kdiff_vs_local_CpG_density(pkl_fp, split=True)
# -*- coding:utf-8 -*-
import math,os,collections
import pandas as pd
import multiprocessing as mp
import matplotlib.pyplot as plt
from localreg import *
import seaborn as sns

data_dir = "../data"
EACH_SUB_FIG_SIZE = 5
D_MAX = 1000
is_inter_with_other_cpg = True
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

def _workers_count():
    try:
        cpu_count = len(os.sched_getaffinity(0))
    except AttributeError:
        cpu_count = os.cpu_count()
    return cpu_count

def read_bed_file_and_store_pos_to_a_struct(bed_tsv_file_path, K_RD, only_within=False, with_dnase_and_nucleo=False):
    """
    Read tsv format bed file (No chr prefix at 1st colum), and store the chr, pos, k or methylation level info into dict
    :param bed_tsv_file_path: tsv format bed file input path
    :return: the dict with chr, pos, k or methylation level info
    """
    df = pd.read_csv(bed_tsv_file_path, sep='\t', header=None).values
    chrs = df[:, 0]
    unique_chr, _ = np.unique(chrs, return_counts=True)
    dict_to_store = {chr_i: {} for chr_i in unique_chr}
    dict_to_store_region_idx = {chr_i: {} for chr_i in unique_chr}
    dict_to_store_dnase_and_nucleo_val = {chr_i: {} for chr_i in unique_chr}
    val_col = 4 if K_RD else 3
    for item in df:
        try:
            if only_within and not with_dnase_and_nucleo:
                chr_i, pos, val, cgi_idx = item[0], item[1], item[val_col], item[-1]
                dict_to_store_region_idx[chr_i][int(pos)] = cgi_idx
            elif with_dnase_and_nucleo:
                chr_i, pos, val, cgi_idx, dnase_val, nucleo_val = item[0], item[1], item[val_col], item[-3], item[-2], item[-1]
                dict_to_store_region_idx[chr_i][int(pos)] = cgi_idx
                dict_to_store_dnase_and_nucleo_val[chr_i][int(pos)] = [float(dnase_val), float(nucleo_val)]
            else:
                chr_i, pos, val = item[0], item[1], item[val_col]
        except Exception as e:
            pass
        dict_to_store[chr_i][int(pos)] = val
    if only_within == False:
        if with_dnase_and_nucleo:
            return [dict_to_store, dict_to_store_dnase_and_nucleo_val]
        else:
            return [dict_to_store]
    else:
        if with_dnase_and_nucleo:
            return [dict_to_store, dict_to_store_region_idx, dict_to_store_dnase_and_nucleo_val]
        else:
            return [dict_to_store, dict_to_store_region_idx]

def filter_d_length_to_generate_CpG_pairs(Chr_CpG_pos_and_methy_list, D_MAX, only_within=False, with_dnase_and_nucleo=False):
    """
    Extract the CpG pairs with distance of d, and there can be CpGs between CpG pairs with d distance
    :param Chr_CpG_pos_and_methy_dict: The dict with chr, pos, k or methylation level info provided
    :param d: the distance that required for CpG pairs
    :return: An array that stores all CpG pairs with distance of d
    """

    CpG_Pairs_Dict = {}
    DNase_Pairs_Dict = {}
    Nucleo_Pairs_Dict = {}
    D_RANGE = range(2, D_MAX)

    Chr_CpG_pos_and_methy_dict = Chr_CpG_pos_and_methy_list[0]
    Chr_CpG_pos_and_cgi_index_dict = {}
    if only_within:
        Chr_CpG_pos_and_cgi_index_dict = Chr_CpG_pos_and_methy_list[1]
        if with_dnase_and_nucleo:
            dict_to_store_dnase_and_nucleo_val = Chr_CpG_pos_and_methy_list[2]
    else:
        if with_dnase_and_nucleo:
            dict_to_store_dnase_and_nucleo_val = Chr_CpG_pos_and_methy_list[1]
    for d in D_RANGE:
        if d % 10 == 0:
            print("Calculating CpG pairs with %d distance" % d)
        for chr_i in Chr_CpG_pos_and_methy_dict.keys():
            CpG_pos_and_methy_dict = Chr_CpG_pos_and_methy_dict[chr_i]
            for key in CpG_pos_and_methy_dict.keys():
                pos = key
                methy_level_1 = CpG_pos_and_methy_dict[key]
                pos_off_by_d = pos + d
                if pos_off_by_d in CpG_pos_and_methy_dict.keys():
                    if (not only_within) or (only_within and (Chr_CpG_pos_and_cgi_index_dict[chr_i][pos_off_by_d] == Chr_CpG_pos_and_cgi_index_dict[chr_i][key])):
                        methy_level_2 = CpG_pos_and_methy_dict[pos_off_by_d]
                        methy_levels_of_CpG_pairs =[methy_level_1, methy_level_2]

                        if with_dnase_and_nucleo and pos:
                            try:
                                p1 = dict_to_store_dnase_and_nucleo_val[chr_i][pos]
                                p2 = dict_to_store_dnase_and_nucleo_val[chr_i][pos_off_by_d]
                                dnase_vals = [p1[0], p2[0]]
                                nucleo_vals = [p1[1], p2[1]]
                                if d not in DNase_Pairs_Dict.keys():
                                    DNase_Pairs_Dict[d]= []
                                    Nucleo_Pairs_Dict[d] =[]
                                DNase_Pairs_Dict[d].append(dnase_vals)
                                Nucleo_Pairs_Dict[d].append(nucleo_vals)
                            except Exception as e:
                                pass
                        if d not in CpG_Pairs_Dict.keys():
                            CpG_Pairs_Dict[d] = []
                        CpG_Pairs_Dict[d].append(methy_levels_of_CpG_pairs)
    if with_dnase_and_nucleo:
        return [CpG_Pairs_Dict, DNase_Pairs_Dict, Nucleo_Pairs_Dict]
    else:
        return [CpG_Pairs_Dict]

def prepare_config_file(CONFIG_FP):
    regions = ["ChromHMM", "Genomic_Regions", "Histone_Modification", "TFBS"]
    ltws = []
    for K in [1]:
        for REGION in regions:
            intersection_dir = os.path.join("../data/K_region_intersect", REGION)
            out_rd_corr_dir = os.path.join(intersection_dir, "K_Rd") if K == 1 else os.path.join(intersection_dir, "Methy_Rd")
            mkdir(out_rd_corr_dir)
            for file_name in os.listdir(intersection_dir):
                if file_name.endswith(".bed"):
                    bed_tsv_file_path = os.path.join(intersection_dir, file_name)
                    out_R_d_correlation_path = os.path.join(out_rd_corr_dir, file_name)
                    ltws.append('\t'.join([bed_tsv_file_path, out_R_d_correlation_path, str(K)]))
    with open(CONFIG_FP, 'w') as CONFIG_file:
        ltw = '\n'.join(ltws) + '\n'
        CONFIG_file.write(ltw)

def filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(Chr_CpG_pos_and_methy_dict, D_MAX):
    """
    Extract the CpG pairs with distance of d, and required no CpGs between CpG pairs with d distance
    :param d: the distance that required for CpG pairs
    :param length_of_od:
    :param od_keys:
    :param od_vals:
    :return:
    """
    CpG_Pairs_list = {}
    Chr_CpG_pos_and_methy_dict=Chr_CpG_pos_and_methy_dict[0]
    for chr_i in Chr_CpG_pos_and_methy_dict.keys():
        CpG_pos_and_methy_dict = Chr_CpG_pos_and_methy_dict[chr_i]
        sorted_struct = collections.OrderedDict(sorted(CpG_pos_and_methy_dict.items()))
        od_keys = list(sorted_struct.keys())
        od_vals = list(sorted_struct.values())
        for i in range(0, len(sorted_struct) - 1):
            # 后一个位点与前一个位点之间的距离为d
            pre_od_key = od_keys[i]
            post_od_key = od_keys[i + 1]
            dist = post_od_key - pre_od_key
            if 0 < dist <D_MAX:
                methy_levels_of_CpG_pairs = [od_vals[i], od_vals[i + 1]]
                if dist not in CpG_Pairs_list.keys():
                    CpG_Pairs_list[dist] = []
                CpG_Pairs_list[dist].append(methy_levels_of_CpG_pairs)
    return [CpG_Pairs_list]

def calc_C_d_by_pearson_correlation(CpG_Pairs_list, with_dnase_and_nucleo=False):
    """
    calc the Pearson correlation of CpG pairs list
    :param CpG_pairs: The input CpG pairs list
    :return: The dict with {d: correlation}
    """
    RD_dict = {}
    if with_dnase_and_nucleo:
        CpG_Pairs_Dict, DNase_Pairs_Dict, Nucleo_Pairs_Dict = CpG_Pairs_list[0], CpG_Pairs_list[1], CpG_Pairs_list[2]
    else:
        CpG_Pairs_Dict = CpG_Pairs_list[0]

    for DISTANCE, CpG_pairs in CpG_Pairs_Dict.items():
        if len(CpG_pairs) > 5:
                CpG_arr = np.array(CpG_pairs)
                try:
                    r_d = np.corrcoef(CpG_arr[:, 0], CpG_arr[:, 1])[0][1]
                    if with_dnase_and_nucleo:
                        mean_dnase = np.mean(np.array(DNase_Pairs_Dict[DISTANCE]).astype(np.float))
                        mean_nucleo = np.mean(np.array(Nucleo_Pairs_Dict[DISTANCE]).astype(np.float))
                        RD_dict[DISTANCE] = [r_d, mean_dnase, mean_nucleo]
                    else:
                        RD_dict[DISTANCE] = r_d
                except Exception as e:
                    pass
    return [RD_dict] if RD_dict else []

def write_RD_into(RD_dict, out_fp, with_dnase_and_nucleo=False):
    if RD_dict:
        RD_dict = RD_dict[0]
        sorted_struct = collections.OrderedDict(sorted(RD_dict.items()))
        od_keys = list(sorted_struct.keys())
        if with_dnase_and_nucleo:
            ltw = "\n".join(["%d\t%.3f\t%.3f\t%.3f" % (d, RD_dict[d][0], RD_dict[d][1], RD_dict[d][2]) for d in od_keys]) + "\n"
        else:
            ltw = "\n".join(["%d\t%.3f" %(d, RD_dict[d]) for d in od_keys]) + "\n"
        with open(out_fp, 'w') as RD_file:
            RD_file.write(ltw)

def calc_correlation(bed_tsv_file_path, out_R_d_correlation_path,  d_max, is_inter_with_other_cpg, K_RD=1, only_within=False, with_dnase_and_nucleo=False):
    Chr_CpG_pos_and_methy_dict = read_bed_file_and_store_pos_to_a_struct(bed_tsv_file_path, K_RD, only_within=only_within, with_dnase_and_nucleo=with_dnase_and_nucleo)
    if is_inter_with_other_cpg:
        CpG_pairs = filter_d_length_to_generate_CpG_pairs(Chr_CpG_pos_and_methy_dict, d_max, only_within=only_within, with_dnase_and_nucleo=with_dnase_and_nucleo)
    else:
        CpG_pairs = filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(Chr_CpG_pos_and_methy_dict, d_max)
    RD_dict = calc_C_d_by_pearson_correlation(CpG_pairs, with_dnase_and_nucleo=with_dnase_and_nucleo)
    write_RD_into(RD_dict, out_R_d_correlation_path, with_dnase_and_nucleo)

def wrapper_for_correlation_analysis(row):
    bed_tsv_file_path = str(row[0])
    out_R_d_correlation_path = str(row[1])
    K_RD = int(row[2])
    print("%s\t%s\t%d" % (bed_tsv_file_path, out_R_d_correlation_path, K_RD))
    calc_correlation(bed_tsv_file_path, out_R_d_correlation_path, D_MAX, is_inter_with_other_cpg, K_RD=K_RD, only_within=False, with_dnase_and_nucleo=True)

def call_for_correlation_analysis(CONFIG_FP):
    num_workers = _workers_count()
    print('Starting call_for_correlation_analysis with', num_workers, 'workers')
    parametersDF = pd.read_csv(CONFIG_FP, index_col=0, sep='\t', header=None, lineterminator='\n')
    tups = parametersDF.itertuples(name=None)
    pool = mp.Pool(processes=num_workers)
    pool.map_async(wrapper_for_correlation_analysis, tups).get()

def plot_local_regression_and_RD(max_d, fig_format="png"):
    regions = ["Genomic_Regions", "Histone_Modification", "TFBS"]# "ChromHMM",
    K_RD = [1]# , 0
    N_COL = 5
    cm = plt.get_cmap('gist_rainbow')
    FIG_DIR = os.path.join(BASE_DIR, "figures")
    for km in K_RD:
        for REGION in regions:
            RD_DIRNAME= "K_Rd" if km == 1 else "Methy_Rd"
            fig_dir = os.path.join(FIG_DIR, RD_DIRNAME)
            mkdir(fig_dir)

            out_rd_corr_dir = os.path.join("../data/K_region_intersect", REGION, RD_DIRNAME)
            file_paths = [os.path.join(out_rd_corr_dir, "%s.bed" % file_name) for file_name in FILE_ORDERED_NAMES[REGION]]
            file_labels = FILE_LABELS[REGION]
            N_FILES = len(file_labels)
            N_ROW = int(math.ceil((N_FILES) / N_COL))
            DNase_Nucleos = ["DNase", "Nuleos_Occupancy"]
            vmins = [0, 120]
            vmaxs = [0.15, 220]
            for dnase_nucleo in range(2):
                fig_fp = os.path.join(fig_dir,  "%s_%s.%s" % (REGION, DNase_Nucleos[dnase_nucleo], fig_format))
                fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
                for j in range(N_FILES):
                    row = j // N_COL
                    col = j % N_COL
                    if N_ROW == 1:
                        ax = axs[col]
                    else:
                        ax = axs[row][col]
                    file_path = file_paths[j]
                    RD_df = pd.read_csv(file_path, sep="\t", header=None).values
                    x = RD_df[:, 0]
                    y = RD_df[:, 1]
                    z = RD_df[:, 2 + dnase_nucleo]
                    try:
                        y2 = localreg(x, y, degree=2, kernel=tricube, width=100)
                        if file_labels[j] == "Genome":
                            sc = ax.scatter(x, y, s=8, c=z, label=file_labels[j], cmap=cm, vmin =vmins[dnase_nucleo], vmax =vmaxs[dnase_nucleo])
                            ax.plot(x, y2, "w-", linewidth=2)
                        else:
                            sc = ax.scatter(x, y, s=8, c=z, label=file_labels[j], cmap=cm, vmin =vmins[dnase_nucleo], vmax =vmaxs[dnase_nucleo])
                            ax.plot(x, y2, "k-", linewidth=2)
                        fig.colorbar(sc, ax=ax)
                    except np.linalg.LinAlgError as e:
                        sns.regplot(x=x, y=y, ax=ax, scatter_kws={'s':8, 'color': cm(1. * j / N_FILES)})#, line_kws ={'color':'black', "lw": 2}
                    ax.set_xticks(range(0, max_d + 1, 200))
                    ax.set_xlim(0, max_d)
                    ax.set_ylim(0, 1.0)
                    ax.set_title(file_labels[j], fontsize=18)
                plt.savefig(fig_fp, dpi=300, bbox_inches='tight', pad_inches=0.1)

if __name__ == "__main__":
    # config_fp = os.path.join(data_dir, "diff_region_with_DNase_and_Nucleo.config")
    # prepare_config_file(config_fp)
    # call_for_correlation_analysis(config_fp)
    plot_local_regression_and_RD(D_MAX)
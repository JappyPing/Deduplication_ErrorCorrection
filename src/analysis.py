# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-03-21 12:00:00
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-06-12 12:15:02

from random import shuffle
from tqdm import tqdm
# from mpire import WorkerPool
from collections import Counter
import pysam
from tqdm import tqdm
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3
import copy
# from matplotlib import pyplot as plt
# from venn import venn
from matplotlib_venn import venn2
from upsetplot import generate_counts, plot
from upsetplot import from_contents
from upsetplot import UpSet
import os
import matplotlib.pyplot as plt
from upsetplot import UpSet
import matplotlib as mpl
import gzip
import numpy as np
from Bio import SeqIO
import math
import seaborn as sns
import matplotlib.ticker as ticker

def parse_file_type(data_set):
    items = data_set.split(".")
    ext = items[-1]
    if ext == 'fa' or ext == 'fasta':
        f_type = 'fasta'
    elif ext == 'fq' or ext == 'fastq':
        f_type = 'fastq'
    elif ext == 'gz':
        f_type = items[-2] + '.' + ext
    return f_type
    
def parse_data(data_set):
    file_type = parse_file_type(data_set)
    if file_type == 'fastq.gz' or file_type == 'fq.gz' or file_type == 'fa.gz' or file_type == 'fasta.gz':
        ff_type = file_type.split('.')[0]
        handle = gzip.open(data_set, 'rt')
        record_iterator = SeqIO.parse(handle, ff_type)
        return record_iterator, ff_type
    else:
        record_iterator = SeqIO.parse(data_set, file_type) 
        return record_iterator, file_type


def faq2seq(f_in):
    seq_lst = []
    # Open the FASTQ file for reading
    with pysam.FastxFile(f_in) as fq:
        # Loop over each record in the file
        for record in fq:
            seq_lst.append(record.sequence)
            # print("fastq")
            # print(record.sequence)
            # Extract the ID and description strings
            # id_parts = record.name.split(" ")
            # id_str = id_parts[0]
            # desc_str = " ".join(id_parts[1:])

            # # Do something with the ID and description
            # id_lst.append(id_str.split(':UMI_')[0])
    return set(seq_lst)   

def bam2seq(f_in):
    seq_lst = []
    # Open the BAM/CRAM file
    bamfile = pysam.AlignmentFile(f_in, "rb")

    # Get all the read sequences in the file
    sequences = []
    for read in bamfile.fetch():
        sequences.append(read.query_sequence)
        # print("bam")
        # print(read.query_sequence)

    # Close the file
    bamfile.close()
    return set(sequences)

def single_end_deduplication(data_base_dir):
    umi_methods = ["AmpUMI", "UMI-tools", "UMIc"] # "UMIc"
    pe_methods = ["Calib", "FastUniq", "Gencore"]
    bam_methods = ["Gencore", "UMI-tools", "umitools"]
    non_umi_methods = ["CD-HIT-DUP", "FastUniq", "pRESTO", "ParDRe", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "Minirmd"]

    dedup_df = pd.DataFrame(columns=["Original", "UMI-tools", "AmpUMI", "UMIc", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "fastp", "pRESTO", "CD-HIT-DUP", "ParDRe", "Minirmd"], dtype=int) # "umitools", 
    new_seq_dedup_df = pd.DataFrame(columns=["Original", "UMI-tools", "AmpUMI", "UMIc", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "fastp", "pRESTO", "CD-HIT-DUP", "ParDRe", "Minirmd"], dtype=int) 
    #######################################
    datasets = ["SRR1543964", "SRR1543965", "SRR1543966","SRR1543967","SRR1543968","SRR1543969","SRR1543970","SRR1543971"]
    # datasets = ["SRR1543964"]
    methods = ["Original", "UMI-tools", "AmpUMI", "UMIc", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "fastp", "pRESTO", "CD-HIT-DUP", "ParDRe", "Minirmd"] #"umitools", 
    mis_methods = ["pRESTO", "CD-HIT-DUP", "Minirmd", "ParDRe"]
    mis_ids = ["0", "1", "2", "3"]

    # Initialize the progress bar
    pbar = tqdm(total=len(datasets)*len(methods))
    for idx in range(len(datasets)):
        dataset = datasets[idx]
        row_lst = []

        new_seq_lst = []

        method2seqs = {}

        ori_path = "/projects/BIOinfo/Jappy/short_read_review/data/hiv_control/fastp_no_umi/"
        ori_data_full_path = ori_path + dataset + ".fastq"
        ori_seqs_set = faq2seq(ori_data_full_path)

        for idx2 in range(len(methods)):
            method = methods[idx2]
            
            if method == "Original":
                ori_path = "/projects/BIOinfo/Jappy/short_read_review/data/hiv_control/fastp_no_umi/"
                data_full_path = ori_path + dataset + ".fastq"
                seqs_lst = faq2seq(data_full_path)
                method2seqs[method] = seqs_lst
                uni_seq_len = len(seqs_lst)
                row_lst.append(uni_seq_len)

                new_len = len(seqs_lst - ori_seqs_set)
                new_seq_lst.append(new_len)
                
            elif method in bam_methods:
                if method == "UMI-tools":
                    data_full_path = data_base_dir + method + "/" + dataset + ".dedup.bam"
                # elif method == "umitools":
                #     data_full_path = data_base_dir + method + "/sorted." + dataset + ".alignment.deumi.sorted.bam"
                # Create the index file
                pysam.index(data_full_path)
                seqs_lst = bam2seq(data_full_path)
                method2seqs[method] = seqs_lst
                uni_seq_len = len(seqs_lst)
                row_lst.append(uni_seq_len)

                new_len = len(seqs_lst - ori_seqs_set)
                new_seq_lst.append(new_len)

            else:
                if method in mis_methods:
                    uni_set = []
                    new_seq_set = []
                    for mis_id in mis_ids:
                        if method == "pRESTO":
                            data_full_path = data_base_dir + method + "/" + mis_id + "/" + dataset + "_collapse-unique.fastq"
                        else:
                            data_full_path = data_base_dir + method + "/" + mis_id + "/" + dataset + ".fastq"                            

                        seqs_lst = faq2seq(data_full_path)
                        method2seqs[method + mis_id] = seqs_lst

                        uni_seq_len = len(seqs_lst)
                        uni_set.append(uni_seq_len)

                        new_len = len(seqs_lst - ori_seqs_set)
                        new_seq_set.append(new_len)
                        
                    new_seq_lst.append(';'.join([str(elem) for elem in new_seq_set]))
                    row_lst.append(';'.join([str(elem) for elem in uni_set]))

                else:
                    if method == "UMIc":
                        data_full_path = data_base_dir + method + "/umi_" + dataset + "_corrected.fastq"
                    elif method == "NGSReadsTreatment":
                        data_full_path = data_base_dir + method + "/" + dataset + "_1_trated.fastq"       
                    elif method == "BioSeqZip":
                        data_full_path = data_base_dir + method + "/" + dataset + ".fq"                                                                  
                    else:
                        data_full_path = data_base_dir + method + "/" + dataset + ".fastq"
                    seqs_lst = faq2seq(data_full_path)
                    method2seqs[method] = seqs_lst                    
                    uni_seq_len = len(seqs_lst)
                    row_lst.append(uni_seq_len)
                    
                    new_len = len(seqs_lst - ori_seqs_set)
                    new_seq_lst.append(new_len)
                    
            pbar.update(1)

        dedup_df.loc[len(dedup_df)] = row_lst
        new_seq_dedup_df.loc[len(new_seq_dedup_df)] = new_seq_lst
        venn_path = "../results/figures/umi_venn/"
        if not os.path.exists(venn_path):
            os.makedirs(venn_path)
        plot_venn3(method2seqs["UMI-tools"], method2seqs["UMIc"], method2seqs["AmpUMI"], ["UMI-tools", "UMIc", "AmpUMI"], venn_path + dataset)
        # using umi-based deduplication results to evaluate non-umi-based methods
        # umi_based_results = {key: method2seqs[key] for key in ["UMI-tools", "UMIc", "AmpUMI"] if key in method2seqs}
        umi_based_results = {}
        umi_based_results["AmpUMI"] = copy.deepcopy(method2seqs["AmpUMI"])
        umi_based_results["UMIc"] = copy.deepcopy(method2seqs["UMIc"])
        umi_based_results["UMI-tools"] = copy.deepcopy(method2seqs["UMI-tools"])
        non_umi_methods = ["NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "fastp", "pRESTO0", "pRESTO1", "pRESTO2", "pRESTO3","CD-HIT-DUP0", "ParDRe0", "Minirmd0", "CD-HIT-DUP1", "ParDRe1", "Minirmd1", "CD-HIT-DUP2", "ParDRe2", "Minirmd2", "CD-HIT-DUP3", "ParDRe3", "Minirmd3"] #"umitools", 
        AmpUMI_y = []
        UMIc_y = []
        UMItools_y = []
        x_stick = []
        for method in non_umi_methods:
            upset_data = copy.deepcopy(umi_based_results)
            del upset_data["UMI-tools"]
            upset_data[method] = copy.deepcopy(method2seqs[method])
            upset_fig_path = "../results/figures/umiEvaNon_umi/" + dataset + "/upsetplot/"
            if not os.path.exists(upset_fig_path):
                os.makedirs(upset_fig_path)            
            plot_upset_plot(upset_data, upset_fig_path+ method + ".png", method, True)
            
            method_unique = method2seqs[method]
            x_stick.append(method)
            AmpUMI_y.append(len(method_unique & method2seqs["AmpUMI"]))
            UMIc_y.append(len(method_unique & method2seqs["UMIc"]))
            UMItools_y.append(len(method_unique & method2seqs["UMI-tools"]))          
        line_graph_path = "../results/figures/umiEvaNon_umi/" + dataset + "/" 
        if not os.path.exists(line_graph_path):
            os.makedirs(line_graph_path)
        plot_heat_map(non_umi_methods, method2seqs, "../results/figures/umiEvaNon_umi/" + dataset + "/heatmap.png", 5, 'tab20c', 90)
        line_chart(x_stick, AmpUMI_y, UMIc_y, UMItools_y, line_graph_path + "line_chart.png")
        # new_ampumi_y = [x / len(umi_based_results["AmpUMI"]) for x in AmpUMI_y]
        # new_umic_y = [x / len(umi_based_results["UMIc"]) for x in UMIc_y]
        # new_umitools_y = [x / len(umi_based_results["UMI-tools"]) for x in UMItools_y]
        
        # stacked_bar_line_chart(x_stick, AmpUMI_y, new_ampumi_y, line_graph_path + "AmpUMI_stacked_bar_line_chart.png")
        # stacked_bar_line_chart(x_stick, UMIc_y, new_umic_y, line_graph_path + "UMIc_stacked_bar_line_chart.png")
        # stacked_bar_line_chart(x_stick, AmpUMI_y, new_umitools_y, line_graph_path + "UMItools_stacked_bar_line_chart.png")
        
    print(dedup_df)
    dedup_df.to_csv("singleEnd_deduplication.csv", index=False, header=True)
    new_seq_dedup_df.to_csv("single_end_new_seq.csv", index=False, header=True)

def plot_heat_map(non_umi_methods, method2seqs, out_fig, fig_height, color_map, rotation):
    data_lst = []
    umic_seqs = method2seqs["UMIc"]
    AmpUMI_seqs = method2seqs["AmpUMI"]
    data_df = pd.DataFrame(index=non_umi_methods, columns=["A", "B", "C", "D", "E", "F", "G"])
    data_df.index = data_df.index.astype(str)
    data_df.columns = data_df.columns.astype(str)
    
    for method in non_umi_methods:
        cur_seqs = method2seqs[method]
        cur_lst = []
        cur_lst.append(len(cur_seqs - umic_seqs - AmpUMI_seqs))
        cur_lst.append(len(AmpUMI_seqs - umic_seqs - cur_seqs))
        cur_lst.append(len(umic_seqs - AmpUMI_seqs - cur_seqs))
        cur_lst.append(len((AmpUMI_seqs & cur_seqs) - umic_seqs))
        cur_lst.append(len((umic_seqs & cur_seqs) - AmpUMI_seqs))
        cur_lst.append(len((umic_seqs & AmpUMI_seqs) - cur_seqs))
        cur_lst.append(len(umic_seqs & AmpUMI_seqs & cur_seqs))
        data_df.loc[method] = cur_lst
    print(data_df)
    data_df = data_df.apply(pd.to_numeric, errors='coerce')

    sns.set()
    num_cols = len(data_df.columns)

    fig, axes = plt.subplots(nrows=1, ncols=num_cols, figsize=(1.2*num_cols, fig_height))

    for i, column in enumerate(data_df.columns):
        cur_cmap = mpl.colormaps[color_map]
        # cmap = sns.color_palette(colormap_names[i], as_cmap=True)
        cbar_kws = {'orientation': 'vertical', 'use_gridspec': True, 'ticks': []}
        if i == 0:
            ax = sns.heatmap(data_df[column].to_frame(),  cmap=cur_cmap, ax=axes[i], annot=True, fmt='.0f', cbar=True, cbar_kws=cbar_kws, annot_kws={'fontsize': 8}) #{'shrink': 1, 'aspect': 30}
        else:
            ax = sns.heatmap(data_df[column].to_frame(),  cmap=cur_cmap, ax=axes[i], yticklabels=False, annot=True, fmt='.0f', cbar=True, cbar_kws=cbar_kws, annot_kws={'fontsize': 8})#{'shrink': 1, 'aspect': 30} #, annot=True
        
        # ax.set_title(column)
        ax.set_xticklabels([])
        if i == 0:
            # ax.set_xticklabels(ax.get_xticklabels(), horizontalalignment='right')
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        # fig.subplots_adjust(left=0.2, right=0.8)
        # fig.subplots_adjust(hspace=0)
        fig.subplots_adjust(left=0.2)
        # Adjust the aspect ratio of the colorbar
        cbar = ax.collections[0].colorbar
        # Set the minimum and maximum values of the colorbar
        cbar_min = data_df[column].min()
        cbar_max = data_df[column].max()

        # Set the position of the colorbar
        # cbar.ax.set_position([0.85, 0.2, 0.03, 0.6])  # Adjust the values as per your requirement

        # Hide the colorbar tick labels
        cbar.ax.set_yticklabels([])

        # Add custom text annotations for the minimum and maximum values
        if cbar_min == 0:
            cbar.ax.text(0.55, -0.05, f'{cbar_min:.0f}', ha='center', va='center',
                        rotation=rotation, transform=cbar.ax.transAxes, fontsize=8)
            cbar.ax.text(0.55, 1.1, f'{cbar_max:.0f}', ha='center', va='center',
                        rotation=rotation, transform=cbar.ax.transAxes, fontsize=8)   
        else:     
            cbar.ax.text(0.55, -0.1, f'{cbar_min:.0f}', ha='center', va='center',
                        rotation=rotation, transform=cbar.ax.transAxes, fontsize=8)
            cbar.ax.text(0.55, 1.1, f'{cbar_max:.0f}', ha='center', va='center',
                        rotation=rotation, transform=cbar.ax.transAxes, fontsize=8)

    plt.tight_layout()

    plt.savefig(out_fig)
    plt.close()

def evaluate_error_correction(raw_dir, cor_dir, dedup_data_base_dir):
    methods = ["Bcool", "Care", "Coral", "Fiona", "Lighter", "Pollux", "RACER", "BFC"] #"Karect", 
    # methods = ["Karect", "Musket", "Coral", "Fiona", "Lighter", "Pollux", "RACER", "BFC"]
    # methods = ["care-cpu", "coral", "fiona", "lighter", "pollux", "bfc"]
    filenames = os.listdir(raw_dir)
    result_df = pd.DataFrame(columns=["Dataset", "Method", "Original Unique Reads Number", "Correct Unique Reads Number", "Corrected Reads Number", "Total Reads Number", "Wrongly Introduced New Reads Number"])

    # filenames = ["SRR1543964.fastq"]
    umi_methods = ["AmpUMI", "UMI-tools", "UMIc"]
    for f_name in filenames:

        umi_dedup_data = {}
        upset_data = {}
        for umi_dedup_method in umi_methods:
            f_name_base = f_name.split(".fastq")[0]
            if umi_dedup_method == "UMIc":
                
                data_full_path = dedup_data_base_dir + umi_dedup_method + "/umi_" + f_name_base + "_corrected.fastq"
                seqs_lst = faq2seq(data_full_path)
                umi_dedup_data[umi_dedup_method] = seqs_lst
                upset_data[umi_dedup_method] = seqs_lst
            elif umi_dedup_method == "UMI-tools":
                data_full_path = dedup_data_base_dir + umi_dedup_method + "/" + f_name_base + ".dedup.bam"        
                pysam.index(data_full_path)
                seqs_lst = bam2seq(data_full_path)
                umi_dedup_data[umi_dedup_method] = seqs_lst
            elif umi_dedup_method == "AmpUMI":
                data_full_path = dedup_data_base_dir + umi_dedup_method + "/" + f_name
                seqs_lst = faq2seq(data_full_path)
                umi_dedup_data[umi_dedup_method] = seqs_lst
                upset_data[umi_dedup_method] = seqs_lst

        ori_path = os.path.join(raw_dir, f_name)
        errCor_upset_data = {}
        for method in methods:
            if method == "Bcool":
                cor_path = os.path.join(cor_dir + method + "/" +f_name, "reads_corrected.fa")
            elif method == "Karect":    
                cor_path = os.path.join(cor_dir + method, "karect_" + f_name)
            elif method == "Lighter":    
                cor_path = os.path.join(cor_dir + method, f_name.split(".")[0] + ".cor.fq")
            else:
                cor_path = os.path.join(cor_dir + method, f_name)
            # result_dir = "./result/" + method + "/"
            # os.system("noise2read -m evaluation -i %s -r %s -d %s" % (ori_path, cor_path, result_dir)) 
            # result_lst = noise2read_entropy(ori_path, cor_path, 60, result_dir, f_name)
            result_lst = get_unique_reads(ori_path, cor_path)
            
            new_lst = result_lst[2:]
            new_lst.insert(0, method)
            new_lst.insert(0, f_name)
            # result_df.append(pd.Series(new_lst, index=result_df.columns), ignore_index=True)
            result_df.loc[len(result_df)] = new_lst

            errCor_upset_data[method] = result_lst[1]
            errCor_upset_data["original"] = result_lst[0]
        upset_fig_path = "../results/figures/ErrCorComparison/"
        if not os.path.exists(upset_fig_path):
            os.makedirs(upset_fig_path) 
        plot_upset_plot(errCor_upset_data, upset_fig_path + f_name + ".png", method, False)
        
        # umi-based deduplication to evaluate error correction
        AmpUMI_y = []
        UMIc_y = []
        UMItools_y = []
        x_stick = []
        for method in methods:
            method_unique = errCor_upset_data[method]

            upset_data[method] = method_unique
            upset_fig_path = "../results/figures/ErrCorComparison/" + f_name + "/upsetplot/"
            if not os.path.exists(upset_fig_path):
                os.makedirs(upset_fig_path)            
            plot_upset_plot(upset_data, upset_fig_path+ method + ".png", method, True)
            del upset_data[method]
            
            x_stick.append(method)
            AmpUMI_y.append(len(method_unique & umi_dedup_data["AmpUMI"]))
            UMIc_y.append(len(method_unique & umi_dedup_data["UMIc"]))
            UMItools_y.append(len(method_unique & umi_dedup_data["UMI-tools"]))          
        line_graph_path = "../results/figures/ErrCorComparison/" + f_name + "/" 
        if not os.path.exists(line_graph_path):
            os.makedirs(line_graph_path)

        line_chart2(x_stick, AmpUMI_y, UMIc_y, UMItools_y, line_graph_path + "line_chart.png")
        merged_data = copy.deepcopy(errCor_upset_data)
        merged_data["AmpUMI"] = umi_dedup_data["AmpUMI"]
        merged_data["UMIc"] = umi_dedup_data["UMIc"]

        plot_heat_map(methods, merged_data, "../results/figures/ErrCorComparison/" + f_name + "/heatmap.png", 2.5, 'coolwarm', 0)
   
    result_df.to_csv("benchmark.csv", index=False)

def get_unique_reads(raw_data, correct_data):
    # read the input using SeqIO
    raw_record_iterator, raw_file_tye = parse_data(raw_data)
    correct_record_iterator, correct_file_tye = parse_data(correct_data)
    raw_seqs = []
    correct_seqs = []

    total_reads_num = 0
    
    raw_record_dict = {}
    correct_record_dict = {}
    id_lst = []
    for raw_item, correct_item in zip(raw_record_iterator, correct_record_iterator):
        raw_id = str(raw_item.id)
        raw_seq = str(raw_item.seq)

        cor_id = str(correct_item.id)
        cor_seq = str(correct_item.seq)
        
        raw_seqs.append(raw_seq)
        correct_seqs.append(cor_seq)
        
        total_reads_num += 1
        
        raw_record_dict[raw_id] = raw_seq
        correct_record_dict[cor_id] = cor_seq
        
        id_lst.append(raw_id)

    raw_unique_reads = set(raw_seqs)
    correct_unique_reads = set(correct_seqs)

    corrected_reads_num = 0
    for item in tqdm(id_lst):
        ori_read = raw_record_dict[item]
        cor_read = correct_record_dict[item]
        if str(cor_read) != str(ori_read):
            corrected_reads_num += 1
            
    new_reads = correct_unique_reads - raw_unique_reads
    new_reads_num = len(new_reads)
    print("Wrongly introduced {} new reads".format(new_reads_num))

    return [raw_unique_reads, correct_unique_reads, len(raw_unique_reads), len(correct_unique_reads), corrected_reads_num, total_reads_num, new_reads_num]

def line_chart2(x_stick, y_1, y_2, y_3, fig_name):

    fig, ax1 = plt.subplots(figsize=(16,10))
    # Set the grid background with a grey rectangle
    sns.set_theme()
    # mpl.rcParams['axes.linewidth'] = 2  # set the default line width of the axis

    ax1.plot(x_stick, y_1, color='#fbb4ae', linewidth=6, marker='o')
    ax1.plot(x_stick, y_2, color='#b3cde3', linewidth=6, marker='s')
    ax1.plot(x_stick, y_3, color='#ccebc5', linewidth=6, marker='^')
    # Add labels to the data points
    for i in range(len(x_stick)):
        plt.text(x_stick[i], y_1[i], f'{y_1[i]}', ha='right')
        plt.text(x_stick[i], y_2[i], f'{y_2[i]}', ha='right')
        plt.text(x_stick[i], y_3[i], f'{y_3[i]}', ha='right')
        
    # ax1.set_title('Stacked Bar Chart')
    # ax1.spines['left'].set_position(('outward', 8))  # set the position of the left axis
    # ax1.spines['bottom'].set_position(('outward', 8))  # set the position of the bottom axis
    ax1.spines['left'].set_visible(False)  # set the position of the left axis
    ax1.spines['bottom'].set_visible(False)  # set the position of the bottom axis
    ax1.spines['right'].set_visible(False)  # hide the right axis
    ax1.spines['top'].set_visible(False)  # hide the top axis

    ax1.legend(fontsize=22, frameon=False, labels=["AmpUMI", "UMIc", "UMI-tools"])
    #  ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), fontsize=18, frameon=False, labels=["AmpUMI", "UMIc", "UMI-tools"])
    plt.xticks(fontsize=22, rotation=90)  # set the font size of the x axis tick labels to 14
    plt.yticks(fontsize=22)  # set the font size of the y axis tick labels to 14
    plt.xlabel('Methods', fontsize=28, fontname="serif")
    ax1.set_ylabel('Common Reads Number', fontsize=28, fontname="serif")
    ax1.spines['left'].set_linewidth(2)
    ax1.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')

    # adjust the layout of the subplots
    # plt.subplots_adjust(wspace=0.3)
    plt.subplots_adjust(bottom=0.2)
    fig.tight_layout()
    plt.savefig(fig_name)
    plt.close(fig)  

def line_chart(x_stick, y_1, y_2, y_3, fig_name):

    fig, ax1 = plt.subplots(figsize=(28,20))

    ax1.plot(x_stick, y_1, color='#fbb4ae', linewidth=4, marker='o')
    ax1.plot(x_stick, y_2, color='#b3cde3', linewidth=4, marker='s')
    ax1.plot(x_stick, y_3, color='#ccebc5', linewidth=4, marker='^')
    # Add labels to the data points
    for i in range(len(x_stick)):
        plt.text(x_stick[i], y_1[i], f'{y_1[i]}', ha='center', va='center', fontsize=22, rotation=45)
        plt.text(x_stick[i], y_2[i], f'{y_2[i]}', ha='center', va='center', fontsize=22, rotation=45)
        plt.text(x_stick[i], y_3[i], f'{y_3[i]}', ha='center', va='center', fontsize=22, rotation=45)
        
    # ax1.set_title('Stacked Bar Chart')
    ax1.spines['left'].set_position(('outward', 8))  # set the position of the left axis
    ax1.spines['bottom'].set_position(('outward', 8))  # set the position of the bottom axis
    # ax1.spines['left'].set_visible(False)  # set the position of the left axis
    # ax1.spines['bottom'].set_visible(False)  # set the position of the bottom axis
    ax1.spines['right'].set_visible(False)  # hide the right axis
    ax1.spines['top'].set_visible(False)  # hide the top axis

    ax1.legend(fontsize=22, frameon=False, labels=["AmpUMI", "UMIc", "UMI-tools"])
    #  ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), fontsize=18, frameon=False, labels=["AmpUMI", "UMIc", "UMI-tools"])
    plt.xticks(fontsize=32, rotation=90)  # set the font size of the x axis tick labels to 14
    plt.yticks(fontsize=32)  # set the font size of the y axis tick labels to 14
    plt.xlabel('Methods', fontsize=32, fontname="serif")
    ax1.set_ylabel('Common Reads Number', fontsize=32, fontname="serif")
    ax1.spines['left'].set_linewidth(2)
    ax1.tick_params(axis='both', labelsize=32)
    ax1.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')

    # adjust the layout of the subplots
    # plt.subplots_adjust(wspace=0.3)
    plt.subplots_adjust(bottom=0.4)
    # plt.grid(True)
    fig.tight_layout()
    plt.savefig(fig_name)
    plt.close(fig)  

def plot_venn3(set1, set2, set3, lab_name_lst, venn_name):
    # Calculate the common sequences
    common_12 = len(set1 & set2)
    common_13 = len(set1 & set3)
    common_23 = len(set2 & set3)
    common_123 = len(set1 & set2 & set3)
    # Set the font name and size for the plot
    plt.rc('font', family='serif', size=14)
    # Create the Venn diagram with common sequence counts
    venn3(subsets=(len(set1)-common_12-common_13-common_123, 
                len(set2)-common_12-common_23-common_123, 
                common_12, 
                len(set3)-common_13-common_23-common_123, 
                common_13, 
                common_23, 
                common_123),
        set_colors=('#1b9e77', '#d95f02', '#7570b3'),
        set_labels=(lab_name_lst[0], lab_name_lst[1], lab_name_lst[2]))
    # Show the plot
    plt.tight_layout()
    plt.savefig(venn_name + ".png")
    plt.close()
    return
#################################################################################

def pair_end_deduplication(data_base_dir):
    umi_methods = ["AmpUMI", "BioSeqZip", "fastp", "UMI-tools"] # "UMIc"
    pe_methods = ["Calib", "FastUniq", "Gencore"]
    bam_methods = ["Gencore", "UMI-tools", "umitools"]
    non_umi_methods = ["CD-HIT-DUP", "FastUniq", "pRESTO", "ParDRe", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "Minirmd"]

    dedup_df = pd.DataFrame(columns=["Original", "AmpUMI", "Calib", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "fastp", "FastUniq", "pRESTO", "CD-HIT-DUP", "ParDRe", "Minirmd"], dtype=int) # "umitools", 
    new_seq_dedup_df = pd.DataFrame(columns=["Original", "AmpUMI", "Calib", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "fastp", "FastUniq", "pRESTO", "CD-HIT-DUP", "ParDRe", "Minirmd"], dtype=int) 
    #######################################
    datasets = ["SRR11207257_1", "SRR11207257_2"]
    methods = ["Original", "AmpUMI", "Calib", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "fastp", "FastUniq", "pRESTO", "CD-HIT-DUP", "ParDRe", "Minirmd"] #"umitools", 

    mis_methods = ["pRESTO", "CD-HIT-DUP", "Minirmd", "ParDRe"]
    mis_ids = ["0", "1", "2", "3"]

    # Initialize the progress bar
    pbar = tqdm(total=len(datasets)*len(methods))
    
    for idx in range(len(datasets)):
        dataset = datasets[idx]
        row_lst = []
        method2seqs = {}
        new_seq_lst = []
        
        ori_path = "/projects/BIOinfo/Jappy/short_read_review/data/hiv/fastp_no_umi/"
        ori_data_full_path = ori_path + dataset + ".fastq"
        ori_seqs_set = faq2seq(ori_data_full_path)
        
        for idx2 in range(len(methods)):
            method = methods[idx2]
            if method == "Original":
                ori_path = "/projects/BIOinfo/Jappy/short_read_review/data/hiv/fastp_no_umi/"
                data_full_path = ori_path + dataset + ".fastq"
                seqs_lst = faq2seq(data_full_path)
                method2seqs[method] = seqs_lst
                uni_seq_len = len(seqs_lst)
                row_lst.append(uni_seq_len)
                new_len = len(seqs_lst - ori_seqs_set)
                new_seq_lst.append(new_len)
            elif method in bam_methods:
                if method == "UMI-tools":
                    data_full_path = data_base_dir + method + "/" + dataset + ".dedup.bam"
                # elif method == "umitools":
                #     data_full_path = data_base_dir + method + "/sorted." + dataset + ".alignment.deumi.sorted.bam"
                # Create the index file
                pysam.index(data_full_path)
                seqs_lst = bam2seq(data_full_path)
                method2seqs[method] = seqs_lst
                uni_seq_len = len(seqs_lst)
                row_lst.append(uni_seq_len)

                new_len = len(seqs_lst - ori_seqs_set)
                new_seq_lst.append(new_len)
                
            else:
                if method in mis_methods:
                    uni_set = []
                    new_seq_set = []
                    for mis_id in mis_ids:
                        data_full_path = data_base_dir + method + "/" + mis_id + "/" + dataset + ".fastq"

                        seqs_lst = faq2seq(data_full_path)
                        method2seqs[method + '_' + mis_id] = seqs_lst

                        uni_seq_len = len(seqs_lst)
                        uni_set.append(uni_seq_len)

                        new_len = len(seqs_lst - ori_seqs_set)
                        new_seq_set.append(new_len)

                    new_seq_lst.append(';'.join([str(elem) for elem in new_seq_set]))
                    row_lst.append(';'.join([str(elem) for elem in uni_set]))
                else:
                    data_full_path = data_base_dir + method + "/" + dataset + ".fastq"
                    seqs_lst = faq2seq(data_full_path)
                    method2seqs[method] = seqs_lst                    
                    uni_seq_len = len(seqs_lst)
                    row_lst.append(uni_seq_len)

                    new_len = len(seqs_lst - ori_seqs_set)
                    new_seq_lst.append(new_len)
            pbar.update(1)
        dedup_df.loc[len(dedup_df)] = row_lst
        new_seq_dedup_df.loc[len(new_seq_dedup_df)] = new_seq_lst
        
        # plot_venn3(method2seqs["UMI-tools"], method2seqs["UMIc"], method2seqs["AmpUMI"], ["UMI-tools", "UMIc", "AmpUMI"], dataset)
        # venn2(method2seqs["UMIc"], method2seqs["AmpUMI"], ["UMIc", "AmpUMI"], dataset)
        
    # print(dedup_df)
    dedup_df.to_csv("pair_end_deduplication.csv", index=False, header=True)
    new_seq_dedup_df.to_csv("pair_end_new_seq.csv", index=False, header=True)

def pair_end_align_deduplication(data_base_dir):

    dedup_df = pd.DataFrame(columns=["Gencore"], dtype=int) # "umitools", , "UMI-tools"

    #######################################
    datasets = ["SRR11207257"]
    methods = ["Gencore"] #"umitools", , "UMI-tools"

    # Initialize the progress bar
    pbar = tqdm(total=len(datasets)*len(methods))
    
    for idx in range(len(datasets)):
        dataset = datasets[idx]
        row_lst = []
        method2seqs = {}
        for idx2 in range(len(methods)):
            method = methods[idx2]
            if method == "UMI-tools":
                data_full_path = data_base_dir + method + "/" + dataset + ".dedup.bam"
            else:
                data_full_path = data_base_dir + method + "/" + dataset + ".bam"
            pysam.index(data_full_path)
            seqs_lst = bam2seq(data_full_path)
            method2seqs[method] = seqs_lst
            uni_seq_len = len(seqs_lst)
            row_lst.append(uni_seq_len)
            pbar.update(1)
        dedup_df.loc[len(dedup_df)] = row_lst
        
        # plot_venn3(method2seqs["UMI-tools"], method2seqs["UMIc"], method2seqs["AmpUMI"], ["UMI-tools", "UMIc", "AmpUMI"], dataset)
        # venn2(method2seqs["UMI-tools"], method2seqs["Gencore"], dataset)
    # print(dedup_df)
    dedup_df.to_csv("pair_end_align_deduplication.csv", index=False, header=True)

def venn2(set1, set2, lab_name_lst, venn_name):
    # Calculate the common sequences
    common = len(set1 & set2)
    plt.rc('font', family='serif', size=14)
    # Create the Venn diagram
    venn2(subsets=(len(set1)-common, len(set2)-common, common),
        set_colors=('#1b9e77', '#d95f02'),
        set_labels=(lab_name_lst[0], lab_name_lst[1]))

    # Show the plot
    plt.savefig(venn_name + ".png")
    plt.close()

def plot_upset_plot(data, fig_name, method, empty_subsets):
    upset_data = from_contents(data)
    fig = plt.figure(figsize=(12, 6))
    # upset = UpSet(upset_data, subset_size='count', show_counts=True, facecolor="#fbb4ae")#orientation='vertical' , include_empty_subsets=True
    upset = UpSet(upset_data, subset_size='count', show_counts=True, facecolor="#fbb4ae", element_size=50, include_empty_subsets=empty_subsets)#orientation='vertical' , include_empty_subsets=True
    # upset.style_subsets(present=[method], facecolor="#b3cde3")
    upset.plot(fig=fig)

    # Set plot title
    # plt.title('UpSet Plot - Common Elements between A, B, C, and D')

    # Show the plot
    plt.savefig(fig_name)
    plt.close()

if __name__ == "__main__":
    dedup_data_base_dir = "/projects/BIOinfo/Jappy/short_read_review/results/deduplication/"
    # single_end_deduplication(dedup_data_base_dir)
    
    raw_dir = "/projects/BIOinfo/Jappy/short_read_review/data/hiv_control/fastp_no_umi/"
    cor_dir = "/projects/BIOinfo/Jappy/short_read_review/results/error_correction/"
    evaluate_error_correction(raw_dir, cor_dir, dedup_data_base_dir) 
    
    # pair_end_deduplication(data_base_dir)
    # pair_end_align_deduplication(data_base_dir)
    # example = generate_counts()
    # print(example)

# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-03-21 12:00:00
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-06-12 12:15:02

# from collections import Counter
import pysam
# from tqdm import tqdm
import pandas as pd
from matplotlib import pyplot as plt
# from matplotlib_venn import venn3
# import copy
# from matplotlib_venn import venn2
# from upsetplot import from_contents
# from upsetplot import UpSet
import os
# import matplotlib.pyplot as plt
# from upsetplot import UpSet
# import matplotlib as mpl
import gzip
from Bio import SeqIO
import seaborn as sns

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

def get_unique_reads(data):

    record_iterator, correct_file_tye = parse_data(data)
    seqs = []
    for item in record_iterator:
        seq = str(item.seq)
        seqs.append(seq)
    unique_reads = set(seqs)
    return unique_reads

def get_unique_reads2(data):
    # Read the content of the file
    with open(data, "r") as file:
        lines = file.readlines()
    seqs = []
    # Process the lines
    for line in lines:
        line = line.strip()
        if line.startswith(">") or line.startswith("@"):
            pass
        else:
            seqs.append(line)
            # print(line)
            # print("---")
    unique_reads = set(seqs)
    return unique_reads

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


def deduplication_after_correction(input_dir, output_dir):
    # mismatch_allowed_methods = ["CD-HIT-DUP", "ParDRe", "Minirmd"]
    mismatch_allowed_methods = ["ParDRe"]
    # correction_methods = ["Bcool", "Care", "Coral", "Fiona", "Lighter", "Pollux", "RACER", "BFC"]
    # correction_methods = ["Bcool", "Coral", "Fiona", "Lighter", "Pollux"]
    # correction_methods = ["Coral", "Fiona"]
    # correction_methods = ["Lighter", "Pollux"]
    correction_methods = ["Care"]#, "BFC", "RACER", "BFC"
    datasets = ["SRR1543964.fastq", "SRR1543965.fastq", "SRR1543966.fastq","SRR1543967.fastq","SRR1543968.fastq","SRR1543969.fastq","SRR1543970.fastq","SRR1543971.fastq"]
    # mis_ids = ["(0)", "(1)", "(2)", "(3)"]
    mis_ids = [1,2]
    # mis_ids = [3]
    for method in mismatch_allowed_methods: 
        for cor_method in correction_methods:
            for data_set in datasets:
                for mis_val in mis_ids:
                    if cor_method == "Bcool":
                        cur_data_path = os.path.join(input_dir + cor_method + "/" + data_set, "reads_corrected.fa")
                    elif cor_method == "Lighter":    
                        cur_data_path = os.path.join(input_dir + cor_method, data_set.split(".")[0] + ".cor.fq")
                    else:
                        cur_data_path = os.path.join(input_dir + cor_method, data_set)
                    out_root_dir = output_dir + method + "/" + str(mis_val) + "/" + cor_method
                    if not os.path.exists(out_root_dir):
                        os.makedirs(out_root_dir)
                    output = os.path.join(out_root_dir, data_set)
                    if method == "CD-HIT-DUP":
                        cd_data_set = data_set.split("fastq")[0] + "fasta"
                        cd_output = os.path.join(out_root_dir, cd_data_set)
                        os.system("cd-hit-dup -i %s -e %s -o %s" % (cur_data_path, str(mis_val), cd_output))
                    elif method == "Minirmd":
                        # print(cur_data_path)
                        os.system("minirmd -i %s -d %s -o %s -t 52" % (cur_data_path, str(mis_val), output))
                    elif method == "ParDRe":
                        bin_method = "/projects/BIOinfo/Jappy/short_read_review/methods/deduplication/ParDRe-rel2.2.5/ParDRe"
                        os.system("mpirun -n 64 %s -i %s -m %s -o %s" % (bin_method, cur_data_path, str(mis_val), output))


def plot_line_chart():
    mismatch_allowed_methods = ["CD-HIT-DUP", "ParDRe", "Minirmd", "no_mis"]
    # mismatch_allowed_methods = ["ParDRe"]
    correction_methods = ["BFC", "Bcool", "Care", "Coral", "Fiona", "Lighter", "Pollux", "RACER"]#, 
    # correction_methods = ["Care", "RACER", "BFC"]
    datasets = ["SRR1543964.fastq", "SRR1543965.fastq", "SRR1543966.fastq","SRR1543967.fastq","SRR1543968.fastq","SRR1543969.fastq","SRR1543970.fastq","SRR1543971.fastq"]
    mis_ids = [1, 2, 3]
    # mis_ids = [3]
    # umi_dedup_methods = ["UMI-tools", "AmpUMI", "UMIc"]
    umi_dedup_methods = ["UMI-tools", "AmpUMI", "UMIc"]
    umi_data_base_dir = "/projects/BIOinfo/Jappy/short_read_review/results/deduplication/"

    de_error_correction_data_base_dir = "/projects/BIOinfo/Jappy/short_read_review/results/de_error_correction/"

    for data_set in datasets:
        method2seqs = {}
        data_set_name = data_set.split(".fastq")[0]
        #umi-based deduplication
        for umi_method in umi_dedup_methods:
            if umi_method == "UMI-tools":
                data_full_path = umi_data_base_dir + umi_method + "/" + data_set_name + ".dedup.bam"                
                pysam.index(data_full_path)
                seqs_lst = bam2seq(data_full_path)
                method2seqs[umi_method] = seqs_lst
            else:
                if umi_method == "UMIc":
                    data_full_path = umi_data_base_dir + umi_method + "/umi_" + data_set_name + "_corrected.fastq"
                
                else:
                    data_full_path = umi_data_base_dir + umi_method + "/" + data_set_name + ".fastq"
                seqs_lst = faq2seq(data_full_path)
                method2seqs[umi_method] = seqs_lst                    
        de_cor_methods_lst = []

        fig, ax = plt.subplots(figsize=(16,10))
        sns.set_theme()
        x_stick = correction_methods
        for mis_val in mis_ids:
            y_ampumi = {}
            y_umic = {}
            y_umitools = {}      
                  
            for mis_method in mismatch_allowed_methods:
                cur_y_amp_umi = []
                cur_y_umic = []
                cur_y_umi_tools = []
                for cor_method in correction_methods:
                    cur_method = mis_method + "_" + str(mis_val) + "_" + cor_method
                    de_cor_methods_lst.append(cur_method)
                    if mis_method == "CD-HIT-DUP":
                        data_path = de_error_correction_data_base_dir + mis_method + "/" + str(mis_val) + "/" + cor_method + "/" + data_set_name + ".fasta"
                    elif mis_method == "no_mis":
                        # data_path = de_error_correction_data_base_dir + "Minirmd/0/" + cor_method + "/" + data_set
                        data_path = de_error_correction_data_base_dir + "CD-HIT-DUP/0/" + cor_method + "/" + data_set_name + ".fasta"
                    else:
                        data_path = de_error_correction_data_base_dir + mis_method + "/" + str(mis_val) + "/" + cor_method + "/" + data_set
                    # if cor_method == "Bcool":
                    #     unique_seq_set = get_unique_reads2(data_path)
                    # else:
                    #     unique_seq_set = get_unique_reads(data_path)
                    unique_seq_set = get_unique_reads2(data_path)
                    method2seqs[cur_method] = unique_seq_set
                    cur_y_amp_umi.append(len(unique_seq_set & method2seqs["AmpUMI"]))
                    cur_y_umic.append(len(unique_seq_set & method2seqs["UMIc"]))
                    cur_y_umi_tools.append(len(unique_seq_set & method2seqs["UMI-tools"]))
                    del unique_seq_set
                y_ampumi[mis_method] = cur_y_amp_umi
                print(f"{mis_method}:")
                print(cur_y_amp_umi)
                y_umic[mis_method] = cur_y_umic
                y_umitools[mis_method] = cur_y_umi_tools

            # Create subplots
            fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(10, 8))
            # "CD-HIT-DUP", "ParDRe", "Minirmd"
            # Plot lines on each subplot
            axes[0].plot(x_stick, y_ampumi["CD-HIT-DUP"], label='CD-HIT-DUP')
            axes[0].plot(x_stick, y_ampumi["ParDRe"], label='ParDRe')
            axes[0].plot(x_stick, y_ampumi["Minirmd"], label='Minirmd')
            axes[0].plot(x_stick, y_ampumi["no_mis"], label='Mismatch=0', color="grey", linestyle="--")
            axes[0].set_title('AmpUMI')
            # axes[0].legend()

            axes[1].plot(x_stick, y_umic["CD-HIT-DUP"])
            axes[1].plot(x_stick, y_umic["ParDRe"])
            axes[1].plot(x_stick, y_umic["Minirmd"])
            axes[1].plot(x_stick, y_umic["no_mis"], label='Mismatch=0', color="grey", linestyle="--")
            axes[1].set_title('UMIc')
            # axes[1].legend()

            axes[2].plot(x_stick, y_umitools["CD-HIT-DUP"])
            axes[2].plot(x_stick, y_umitools["ParDRe"])
            axes[2].plot(x_stick, y_umitools["Minirmd"])
            axes[2].plot(x_stick, y_umitools["no_mis"], label='Mismatch=0', color="grey", linestyle="--")
            axes[2].set_title('UMI-tools')
            # axes[2].legend()
            # Add a single legend at the upper center with one row
            # fig.legend(loc='upper right', bbox_to_anchor=(0.5, 1.0), ncol=1)#bbox_to_anchor=(0.5, 1.1),
            fig.legend(['CD-HIT-DUP', 'ParDRe', 'Minirmd', "Mismatch=0"], loc='upper right', ncol=2, fontsize='small')
            # Adjust layout to add more spacing between legend and subplot titles
            # plt.subplots_adjust(top=0.8)  # Increase the top margin
            # Add annotations for values
            for ax in axes:
                # print(ax)
                for line in ax.lines:
                    # print(line)
                    # if line.properties().get('label') == 'CD-HIT-DUP' or line.properties().get('label') == '_child0': #
                    if line.properties().get('label') == 'Mismatch=0'  or line.properties().get('label') == '_child3':
                        for i in range(len(x_stick)):
                            ax.annotate(f'{line.get_ydata()[i]:.0f}',
                                        (x_stick[i], line.get_ydata()[i]),
                                        fontsize=8)
                            

            # Set common labels and show the plot
            plt.xlabel('Error Correction Methods')
            plt.tight_layout()
            # plt.show()
            # plt.subplots_adjust(bottom=0.2)
            fig.tight_layout()
            
            line_graph_path = "../results/figures/dedup_error_correction/" + str(mis_val) +"/"
            if not os.path.exists(line_graph_path):
                os.makedirs(line_graph_path)        
            
            plt.savefig(line_graph_path + data_set_name + "_line_chart.png")
            plt.close(fig)  

        #         ax.plot(x_stick, y_amp_umi, color='#fbb4ae', linewidth=2, marker='o')
        #         ax.plot(x_stick, y_umic, color='#b3cde3', linewidth=2, marker='s')
        #         ax.plot(x_stick, y_umi_tools, color='#ccebc5', linewidth=2, marker='^')                
        #         for i in range(len(x_stick)):
        #             plt.text(x_stick[i], y_amp_umi[i], f'{y_amp_umi[i]}', ha='right')
        #             plt.text(x_stick[i], y_umic[i], f'{y_umic[i]}', ha='right')
        #             plt.text(x_stick[i], y_umi_tools[i], f'{y_umi_tools[i]}', ha='right')
        # # umi-based deduplication to evaluate error correction
        # ax.spines['left'].set_visible(False)  # set the position of the left axis
        # ax.spines['bottom'].set_visible(False)  # set the position of the bottom axis
        # ax.spines['right'].set_visible(False)  # hide the right axis
        # ax.spines['top'].set_visible(False)  # hide the top axis
        # plt.xticks(fontsize=22, rotation=90)  # set the font size of the x axis tick labels to 14
        # plt.yticks(fontsize=22)  # set the font size of the y axis tick labels to 14
        # plt.xlabel('Methods', fontsize=28, fontname="serif")
        # ax.set_ylabel('Common Reads Number', fontsize=28, fontname="serif")
        # ax.spines['left'].set_linewidth(2)
        # ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')

        # adjust the layout of the subplots
        # plt.subplots_adjust(bottom=0.2)
        # fig.tight_layout()
        
        # line_graph_path = "../results/figures/dedup_error_correction/" 
        # if not os.path.exists(line_graph_path):
        #     os.makedirs(line_graph_path)        
        
        # plt.savefig(line_graph_path + data_set_name + "_line_chart.png")
        # plt.close(fig)  
    return


if __name__ == "__main__":
    input_dir = "/projects/BIOinfo/Jappy/short_read_review/results/error_correction/"
    output_dir = "/projects/BIOinfo/Jappy/short_read_review/results/de_error_correction/new/"
    # deduplication_after_correction(input_dir, output_dir)
    plot_line_chart()


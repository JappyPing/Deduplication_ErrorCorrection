# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-04-13 09:25:44
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-13 16:17:22
import os
from random import shuffle
from tqdm import tqdm
from mpire import WorkerPool
import numpy as np
from Bio import SeqIO
import gzip
from collections import Counter
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

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

def gain2heatmap(correct2raw_diff, result_dir, prefix):
    
    x_len = len(correct2raw_diff)
    x_len_sqrt = round(math.sqrt(x_len))
    xx = x_len_sqrt * x_len_sqrt
    val = x_len - xx
    if xx == x_len:
        m = x_len_sqrt
        n = x_len_sqrt
        x = np.array((correct2raw_diff))
        x_res=x.reshape(m, n)

    elif val < x_len_sqrt:
        m = x_len_sqrt
        n = x_len_sqrt + 1
        remain = (m - val) * [0]
        x = np.array((correct2raw_diff + remain))
        x_res = x.reshape(m, n)

    elif val > x_len_sqrt:
        m = x_len_sqrt
        n = x_len_sqrt + 2
        remain = (val - m) * [0]
        x = np.array((correct2raw_diff + remain))
        x_res=x.reshape(m, n)

    fig, ax = plt.subplots(figsize=(6, 6))

    cur_cmap = mpl.cm.get_cmap('tab20c').copy()#
    # cur_cmap.set_over('red')
    cur_cmap.set_under('red')
    cur_cmap.set_bad(color='red')
    image = ax.imshow(x_res, cmap=cur_cmap) #interpolation='none', , vmax=thre_max
    # image = ax.pcolor(x_res, cmap=cmap, vmin=thre_min, antialiased=True)
    cbar = fig.colorbar(image, extend='min', shrink=0.8)
    cbar.cmap.set_under('red')
    # cbar.cmap.set_over('red')
    fig.tight_layout()
    fig.savefig(os.path.join(result_dir, prefix + '_information_gain.png'), transparent=True)

    return

def entropy_item(total_freq, freq):
    if freq > 0:
        p = freq / total_freq
        result = - p * math.log2(p)
    return result

def noise2read_entropy(raw_data, correct_data, num_workers, result_dir, prefix):

    # read the input using SeqIO
    raw_record_iterator, raw_file_tye = parse_data(raw_data)
    correct_record_iterator, correct_file_tye = parse_data(correct_data)
    raw_seqs = []
    correct_seqs = []

    # raw_len2seqs_dict = {}
    # raw_len_lst = []
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
        
        # ll = len(raw_seq)
        # raw_len_lst.append(ll)
        # raw_len2seqs_dict.setdefault(ll, []).append(raw_seq)
        total_reads_num += 1
        
        raw_record_dict[raw_id] = raw_seq
        correct_record_dict[cor_id] = cor_seq
        
        id_lst.append(raw_id)

    corrected_reads_num = 0

    for item in tqdm(id_lst):
        ori_read = raw_record_dict[item]
        cor_read = correct_record_dict[item]
        if str(cor_read) != str(ori_read):
            corrected_reads_num += 1
            
    raw_read2count = Counter(raw_seqs)
    correct_read2count = Counter(correct_seqs)

    raw_unique_reads = set(raw_read2count.keys())
    correct_unique_reads = set(correct_read2count.keys())

    frequent_reads = set([k for k, v in raw_read2count.items() if v > 5])

    # raw dataset
    non_frequent_raw_reads = raw_unique_reads - frequent_reads
    raw_entropy_items = []
    for read in non_frequent_raw_reads:
        raw_entropy_items.append(raw_read2count[read])

    raw_nonFre_reads_total_num = sum(raw_entropy_items)
    # raw entropy
    with WorkerPool(num_workers, shared_objects=raw_nonFre_reads_total_num, start_method='fork') as pool:
        raw_entropy_lst = pool.map(entropy_item, raw_entropy_items)
    raw_entropy = sum(raw_entropy_lst) 

    # correct dateset
    non_frequent_correct_reads = correct_unique_reads - frequent_reads
    correct_entropy_items = []
    for read in non_frequent_correct_reads:
        correct_entropy_items.append(correct_read2count[read])

    # correct entropy
    correct_nonFre_reads_total_num = sum(correct_entropy_items)

    with WorkerPool(num_workers, shared_objects=correct_nonFre_reads_total_num, start_method='fork') as pool:
        correct_entropy_lst = pool.map(entropy_item, correct_entropy_items) 
    correct_entropy = sum(correct_entropy_lst)
    ##################################################################################
    #information gain (\delta I) heatmap
    new_reads = correct_unique_reads - raw_unique_reads
    new_reads_num = len(new_reads)
    print("Wrongly introduced {} new reads".format(new_reads_num))

    raw_kept_counts = []
    correct_kept_counts = []
    kept_reads = correct_unique_reads & raw_unique_reads
    for read in kept_reads:
        correct_kept_counts.append(correct_read2count[read])
        raw_kept_counts.append(raw_read2count[read])

    with WorkerPool(num_workers, shared_objects=total_reads_num, start_method='fork') as pool:
        raw_kept_entropy_lst = pool.map(entropy_item, raw_kept_counts)

    with WorkerPool(num_workers, shared_objects=total_reads_num, start_method='fork') as pool:
        correct_kept_entropy_lst = pool.map(entropy_item, correct_kept_counts) 

    raw_removed_reads = raw_unique_reads - correct_unique_reads
    raw_removed_items = []
    for read in raw_removed_reads:
        raw_removed_items.append(raw_read2count[read])

    with WorkerPool(num_workers, shared_objects=total_reads_num, start_method='fork') as pool:
        raw_removed_entropy_items_lst = pool.map(entropy_item, raw_removed_items)
    for i in raw_removed_entropy_items_lst:
        if i <=0:
            print('Warning')
    entropy_item_lst = []
    for i, j in zip(raw_kept_entropy_lst, correct_kept_entropy_lst):
        entropy_item_lst.append(i - j)
    entropy_item_lst.extend(raw_removed_entropy_items_lst)
    if new_reads_num > 0:
        entropy_item_lst.extend([np.nan] * new_reads_num)
    shuffle(entropy_item_lst)
    gain2heatmap(entropy_item_lst, result_dir, prefix)
    # return [raw_entropy, correct_entropy]
    return [raw_unique_reads, correct_unique_reads, len(raw_unique_reads), len(correct_unique_reads), corrected_reads_num, total_reads_num, new_reads_num, raw_entropy, correct_entropy]

def evaluate_error_correction(raw_dir, cor_dir):
    datasets = ["SRR1543964", "SRR1543965", "SRR1543966","SRR1543967","SRR1543968","SRR1543969","SRR1543970","SRR1543971"]
    methods = ["bcool", "bfc", "care-cpu", "coral", "fiona", "karect", "lighter", "musket", "pollux", "RACER"]
    filenames = os.listdir(raw_dir)
    result_df = pd.DataFrame(columns=["Dataset", "Original Unique Reads Number", "Correct Unique Reads Number", "Corrected Reads Number", "Total Reads Number", "Wrongly Introduced New Reads Number", "Raw Entropy", "Correct Entropy"])
    for f_name in filenames:

        ori_path = os.path.join(raw_dir, f_name)
        for method in methods:
            if method == "bcool":
                cor_path = os.path.join(cor_dir + method + "/" +f_name, "reads_corrected.fa")
            elif method == "karect":    
                cor_path = os.path.join(cor_dir + method, "karect_" + f_name)
            else:
                cor_path = os.path.join(cor_dir + method, f_name)
            result_dir = "./result/" + method + "/"
            # os.system("noise2read -m evaluation -i %s -r %s -d %s" % (ori_path, cor_path, result_dir)) 
            result_lst = noise2read_entropy(ori_path, cor_path, 60, result_dir, f_name)
            new_lst = result_lst[2:]
            new_lst.insert(0, f_name)
            # result_df.append(pd.Series(new_lst, index=result_df.columns), ignore_index=True)
            result_df.loc[len(result_df)] = new_lst
    result_df.to_csv("benchmark.csv", index=False)

if __name__ == "__main__":
    raw_dir = "/projects/BIOinfo/Jappy/review/data/hiv_control/fastp_no_umi/"
    cor_dir = "/projects/BIOinfo/Jappy/review/result/correction/"
    evaluate_error_correction(raw_dir, cor_dir)

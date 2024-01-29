# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-03-18 21:27:55
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-06-08 23:26:14

# from evaluation import DataAnalysis
import os
from Bio import SeqIO
import gzip
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm
from mpire import WorkerPool

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
    
def parse_data_index(data_set):
    file_type = parse_file_type(data_set)
    if file_type == 'fastq.gz' or file_type == 'fq.gz' or file_type == 'fa.gz' or file_type == 'fasta.gz':
        ff_type = file_type.split('.')[0]
        handle = gzip.open(data_set, 'rt')
        records = SeqIO.index(handle, ff_type)
        return records, ff_type
    else:
        records = SeqIO.index(data_set, file_type)     
        return records, file_type

def parse_data_dict(data_set):
    file_type = parse_file_type(data_set)
    if file_type == 'fastq.gz' or file_type == 'fq.gz' or file_type == 'fa.gz' or file_type == 'fasta.gz':
        ff_type = file_type.split('.')[0]
        handle = gzip.open(data_set, 'rt')
        records_dict = SeqIO.to_dict(SeqIO.parse(handle, ff_type))
        return records_dict, ff_type
    else:
        records_dict = SeqIO.to_dict(SeqIO.parse(data_set, file_type))   
        return records_dict, file_type

def covert_fq2fa(fin, fout):
    records_dict, _ = parse_data_dict(fin)
    records, _ = parse_data_index(fin)
    name_lst = list(records)
    fa_records = []
    with WorkerPool(30, shared_objects=records_dict, start_method='fork') as pool:
        with tqdm(total=len(name_lst), desc="Convert fq to fa") as pbar:   
            for tmp_rec in pool.imap(fq2fa, name_lst):
                fa_records.append(tmp_rec)
                pbar.update() 
    with open(fout, "w") as handle:
        SeqIO.write(fa_records, handle, 'fasta')

def fq2fa(shared_objects, name):
    return SeqRecord(shared_objects[name].seq, id=shared_objects[name].id, description=shared_objects[name].description)

def run(input_dir, output_dir):
    # raw_dir = base_dir + "raw/"
    # true_dir = base_dir + "true/"
    # output_dir = base_dir + "evaluation/bcool/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_names = os.listdir(input_dir)
    raw_fa_dir = output_dir +  "raw_fa/"
    if not os.path.exists(raw_fa_dir):
        os.mkdir(raw_fa_dir)

    # ###############################################################
    for file_name in file_names:
        raw = os.path.join(input_dir, file_name)
        base_name = file_name.split(".fastq")[0]
        raw_fa = raw_fa_dir + base_name + ".fasta"
        covert_fq2fa(raw, raw_fa)   

        correct_dir = os.path.join(output_dir, file_name + "/input/")
        if not os.path.exists(correct_dir):
            os.makedirs(correct_dir)
        os.system("bcool -t 55 -u %s -o %s" % (raw_fa, correct_dir))



if __name__ == "__main__":
    input_dir = "/projects/BIOinfo/Jappy/review/data/hiv_control/fastp_no_umi/"
    output_dir = "/projects/BIOinfo/Jappy/review/result/correction/bcool/"
    run(input_dir, output_dir)

    # input_dir2 = "/projects/BIOinfo/Jappy/review/data/hiv/fastp_no_umi/"
    # run(input_dir2, output_dir)

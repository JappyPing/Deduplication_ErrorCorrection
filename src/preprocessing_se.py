# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-03-13 13:00:44
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-03-13 21:41:13

# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-03-13 12:19:11
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-03-13 12:55:46
from Bio import SeqIO
import gzip
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# from skbio import io
import pysam
# import argparse

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

def preprocessing_data(original_f, trimmed_f, out_f):
    # ori_rec, f_tye = parse_data(original_f)
    record_dict, ori_file_type = parse_data_dict(original_f)
    # trimmed_rec, f_tye = parse_data(trimmed_f)
    print("test")
    id_lst = []
    # for item in tqdm(trimmed_rec):
    #     id_name = str(item.id).split(':UMI_')[0]  
    #     id_lst.append(id_name)  

    # Open the FASTQ file for reading
    with pysam.FastxFile(trimmed_f) as fq:

        # Loop over each record in the file
        for record in fq:

            # Extract the ID and description strings
            id_parts = record.name.split(" ")
            id_str = id_parts[0]
            desc_str = " ".join(id_parts[1:])

            # Do something with the ID and description
            id_lst.append(id_str.split(':UMI_')[0])

    # print('test')
    new_records = []
    # new_fastp_records = []
    for id_n in tqdm(id_lst):
        read = record_dict[id_n].seq
        new_seq = read[0:12] + read[24:]
        # print(len(new_seq))

        qual = {}
        val = record_dict[id_n].letter_annotations['phred_quality'] 
        qual['phred_quality'] = val[0:12] + val[24:]      
        new_rec = SeqRecord(Seq(new_seq), id=id_n, description=record_dict[id_n].description, letter_annotations=qual)
        # print(len(qual['phred_quality']))  
        new_records.append(new_rec)  

    with open(out_f, "w") as handle:
        SeqIO.write(new_records, handle, ori_file_type)

def preprocessing_no_umi(original_f, trimmed_f, out_f):
    # ori_rec, f_tye = parse_data(original_f)
    record_dict, ori_file_type = parse_data_dict(original_f)
    # trimmed_rec, f_tye = parse_data(trimmed_f)
    print("test")
    id_lst = []
    # for item in tqdm(trimmed_rec):
    #     id_name = str(item.id).split(':UMI_')[0]  
    #     id_lst.append(id_name)  

    # Open the FASTQ file for reading
    with pysam.FastxFile(trimmed_f) as fq:

        # Loop over each record in the file
        for record in fq:

            # Extract the ID and description strings
            id_parts = record.name.split(" ")
            id_str = id_parts[0]
            # Do something with the ID and description
            id_lst.append(id_str.split(':UMI_')[0])

    print('test')
    new_records = []
    # new_fastp_records = []
    for id_n in tqdm(id_lst):
        read = record_dict[id_n].seq
        new_seq = read[24:]
        # print(len(new_seq))

        qual = {}
        val = record_dict[id_n].letter_annotations['phred_quality'] 
        qual['phred_quality'] = val[24:]      
        new_rec = SeqRecord(Seq(new_seq), id=id_n, description=record_dict[id_n].description, letter_annotations=qual)
        # print(len(qual['phred_quality']))  
        new_records.append(new_rec)  

    with open(out_f, "w") as handle:
        SeqIO.write(new_records, handle, ori_file_type)

if __name__ == "__main__":
    # original_f = "./original/SRR1543966.fastq"
    # trimmed_f = "./fastp_umi_in_name/umi_SRR1543966.fastq"
    # out_f = "./fastp_umi_in_read/umi_SRR1543966.fastq"

    # Create an argument parser object
    # parser = argparse.ArgumentParser(description='preprocessing data')

    # # Add an argument to the parser
    # parser.add_argument('-s', type=str, required=True,
    #                     help='Path to input original file')
    # parser.add_argument('-i', type=str, required=True,
    #                     help='Path to input fastp processed file')
    # parser.add_argument('-o', type=str, required=True,
    #                     help='Path to output file')
    # # Parse the command-line arguments
    # args = parser.parse_args()    
    # preprocessing_data(args.s, args.i,args.o)
    name_lst = ["SRR1543964.fastq", "SRR1543965.fastq", "SRR1543966.fastq", "SRR1543967.fastq","SRR1543968.fastq","SRR1543969.fastq", "SRR1543970.fastq", "SRR1543971.fastq"]
    ori_path = "/projects/BIOinfo/Jappy/review/data/hiv_control/original/"
    out_dir = "/projects/BIOinfo/Jappy/review/data/hiv_control/fastp_no_umi/"
    trimmed_dir = "/projects/BIOinfo/Jappy/review/data/hiv_control/fastp_umi_in_name/"
    for name in name_lst:
        preprocessing_no_umi(ori_path+name, trimmed_dir + "umi_"+name, out_dir+name)

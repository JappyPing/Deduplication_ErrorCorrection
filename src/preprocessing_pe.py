# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-03-13 12:19:11
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-03-13 22:00:20
from Bio import SeqIO
import gzip
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pysam
from tqdm import tqdm

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
    
    id_lst = []
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

    new_records = []
    for id_n in tqdm(id_lst):
        new_rec = SeqRecord(Seq(record_dict[id_n].seq), id=id_n, description=record_dict[id_n].description, letter_annotations=record_dict[id_n].letter_annotations)  
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
    original_f = "./original/SRR11207257_1.fastq"
    trimmed_f = "./umi_in_name/trimmed_umi_SRR11207257_1.fastq"
    # out_f = "./umi_in_read/SRR11207257_1.fastq"
    out_f = "./fastp_no_umi/SRR11207257_1.fastq"
    # preprocessing_data(original_f, trimmed_f, out_f)
    preprocessing_no_umi(original_f, trimmed_f, out_f)

    original_f2 = "./original/SRR11207257_2.fastq"
    trimmed_f2 = "./umi_in_name/trimmed_umi_SRR11207257_2.fastq"
    out_f2 = "./fastp_no_umi/SRR11207257_2.fastq"
    # preprocessing_data(original_f2, trimmed_f2, out_f2)

    preprocessing_no_umi(original_f2, trimmed_f2, out_f2)
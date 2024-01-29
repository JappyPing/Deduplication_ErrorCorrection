# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2022-12-19 18:35:42
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-06-11 19:48:04
'''
Author: Pengyao PING
Date: 2022-08-01 17:22:47
LastEditors: Pengyao PING
LastEditTime: 2022-10-19 23:39:55
Email: Pengyao.Ping@student.uts.edu.au
Description: 
'''
import os

# def correction(method, raw_dir, correct_base_dir):
#     file_names = os.listdir(raw_dir)
#     bin_method = "./errorCorrection/" + method
#     correct_dir = correct_base_dir + method + "/"
#     if not os.path.exists(correct_dir):
#         os.makedirs(correct_dir)

#     for file_name in file_names:
#         raw = os.path.join(raw_dir, file_name)
#         correct = os.path.join(correct_dir, file_name)
#         gl = 9600
    
#         if method == "bfc":
#             k = 30
#             os.system("%s -t26 -k %s %s > %s" % (bin_method, k, raw, correct))
#         elif method == "coral":
#             k = 26     
#             os.system("%s -fq %s -p 60 -o %s -k %s" % (bin_method, raw, correct, k))
#         elif method == "fiona":
#             os.system("%s -nt 60 -g %s %s %s" % (bin_method, gl, raw, correct))
#             # os.system("%s -nt 30 %s %s" % (bin_method, raw, correct))
#         elif method == "lighter":
#             k = 30
#             os.system("%s -t 60 -r %s -K %s 9600 -od %s" % (bin_method, raw, k, correct_dir))   
#         elif method == "musket":
#             k = 30
#             os.system("LD_LIBRARY_PATH=/home/pping/anaconda3/envs/accHiFi/lib %s -k %s gl -o %s -p 20 %s" % (bin_method, k, correct,raw))        
#         elif method == "pollux":
#             k = 28
#             os.system("%s -i %s -o %s -k %s" % (bin_method, raw, correct_dir, k))   
#             os.system("cat %s %s > %s" % (correct_dir + file_name + ".corrected", correct_dir + file_name + ".low", correct_dir + file_name))
#         elif method == "RACER":
#             os.system("%s %s %s %s" % (bin_method, raw, correct, gl))   
#         elif method == "karect":
#             os.system("%s -correct -threads=60 -matchtype=edit -celltype=haploid  -minoverlap=120 -kmer=14 -kmererrors=2 -inputfile=%s -resultdir=%s" % (bin_method, raw, correct_dir)) 
#         elif method == "care-cpu":
#             coverage = 64
#             os.system("LD_LIBRARY_PATH=/home/pping/anaconda3/envs/accHiFi/lib %s -p -i %s -c %s -o %s -d %s -m 27G -t 60 --pairmode SE -h 48" % (bin_method, raw, coverage, file_name, correct_dir)) 
#             # os.system("%s -c %s -i %s -o %s -t 60 -d %s --pairmode se" % (bin_method, coverage, raw, correct, correct_dir))             
#     return  


def correction(method, raw_dir, correct_base_dir):
# def correction(method, raw, correct_base_dir):
    file_names = os.listdir(raw_dir)
    # file_name = raw.split("/")[-1]
    bin_method = "./errorCorrection/" + method
    correct_dir = correct_base_dir + method + "/"
    if not os.path.exists(correct_dir):
        os.makedirs(correct_dir)

    for file_name in file_names:

        raw = os.path.join(raw_dir, file_name)
        correct = os.path.join(correct_dir, file_name)
        genome_len = 10000
        if method == "bfc":
            k = 30
            os.system("%s -t26 -k %s %s > %s" % (bin_method, k, raw, correct))
        elif method == "coral":
            if "SRR1543965" in file_name:
                k = 28
            else:
                k = 30
            os.system("%s -fq %s -p 26 -o %s -k %s" % (bin_method, raw, correct, k))
        elif method == "fiona":
            os.system("%s -nt 26 -g %s %s %s" % (bin_method, genome_len, raw, correct))
        elif method == "lighter":
            if "SRR1543967" in file_name or "SRR1543970" in file_name:
                k = 28
            elif "SRR1543968" in file_name:
                k = 26
            else:
                k = 30
            os.system("%s -t 26 -r %s -K %s %s -od %s" % (bin_method, raw, k, genome_len, correct_dir))   
        elif method == "musket":
            k = 20
            os.system("%s -k %s genome_len -o %s -p 20 %s" % (bin_method, k, correct,raw))        
        elif method == "pollux":
            k = 20
            os.system("%s -i %s -o %s -k %s" % (bin_method, raw, correct_dir, k))   
            os.system("cat %s %s > %s" % (correct_dir + file_name + ".corrected", correct_dir + file_name + ".low", correct_dir + file_name))
        elif method == "RACER":
            os.system("%s %s %s %s" % (bin_method, raw, correct, genome_len))   
        elif method == "karect":
            os.system("%s -correct -threads=26 -matchtype=edit -celltype=haploid  -minoverlap=120 -kmer=14 -kmererrors=2 -inputfile=%s -resultdir=%s" % (bin_method, raw, correct_dir)) 
        elif method == "care-cpu":
            coverage = 60
            os.system("%s -p -i %s -c %s -o %s -d %s -m 27G -t 26 --pairmode SE -h 48" % (bin_method, raw, coverage, file_name, correct_dir)) 
    
    return  

if __name__ == "__main__":
    raw_dir = "/projects/BIOinfo/Jappy/review/data/hiv_control/fastp_no_umi/"
    # raw_dir = "/projects/BIOinfo/Jappy/review/data/hiv/fastp_no_umi/"

    correct_base_dir = "/projects/BIOinfo/Jappy/review/result/correction/"
    
    # methods = ["bfc", "coral", "musket", "care-cpu", "musket", "fiona", "karect", "lighter", "pollux"]
    methods = ["karect"]
    # file_name = "SRR1543964.fastq"
    # file_name = "SRR1543965.fastq"
    # file_name = "SRR1543966.fastq"
    # file_name = "SRR1543967.fastq"
    # file_name = "SRR1543968.fastq"
    # file_name = "SRR1543969.fastq"
    # file_name = "SRR1543970.fastq"
    # file_name = "SRR1543971.fastq"
    # raw = raw_dir + file_name
    for method in methods:
        correction(method, raw_dir, correct_base_dir)


#!/bin/bash
# @Author: Pengyao Ping
# @Date:   2023-03-15 10:02:53
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-03-24 17:14:36
# mismathes="0 1 2 3"

path="/projects/BIOinfo/Jappy/review/data/hiv/fastp_no_umi/"
out_path_base="/projects/BIOinfo/Jappy/review/result/deduplication/ParDRe/"

data_set1_fullpath=$(readlink -f "${path}SRR11207257_1.fastq")
data_set2_fullpath=$(readlink -f "${path}SRR11207257_2.fastq")
echo ${data_set_fullpath}

out_path="${out_path_base}mis2/"
output_file1_fullpath=$(readlink -f "${out_path}SRR11207257_1.fastq")
output_file2_fullpath=$(readlink -f "${out_path}SRR11207257_2.fastq")
mpirun -n 16 ./ParDRe -i ${data_set1_fullpath} -p ${data_set2_fullpath} -m 2 -o ${output_file1_fullpath} -r ${output_file2_fullpath}


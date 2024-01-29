#!/bin/bash
# @Author: Pengyao Ping
# @Date:   2023-03-14 10:37:40
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-03-14 11:11:01

datasets="SRR1543964.fastq SRR1543965.fastq SRR1543966.fastq SRR1543967.fastq SRR1543968.fastq SRR1543969.fastq SRR1543970.fastq SRR1543971.fastq"

mismathes="0 1 2 3"

path="/projects/BIOinfo/Jappy/review/data/hiv_control/fastp_no_umi/"
out_path_base="/projects/BIOinfo/Jappy/review/result/deduplication/ParDRe/"

for data_set in $datasets
do
    
    srr=${data_set%.fastq}
    
    data_set_fullpath=$(readlink -f "${path}${data_set}")
    echo ${data_set_fullpath}

    for mis_n in $mismathes
    do
        out_path="${out_path_base}mis${mis_n}/"
        output_file_fullpath=$(readlink -f "${out_path}${srr}.fastq")
        mpirun -n 8 ./ParDRe -i ${data_set_fullpath} -m ${mis_n} -o ${output_file_fullpath}
    done
done
#!/bin/bash
# @Author: Pengyao Ping
# @Date:   2023-03-14 10:37:40
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-03-14 11:55:24

# datasets="SRR1543964.fastq SRR1543965.fastq SRR1543966.fastq SRR1543967.fastq SRR1543968.fastq SRR1543969.fastq SRR1543970.fastq SRR1543971.fastq"
datasets="SRR11207257_1.fastq SRR11207257_2.fastq"


# path="/projects/BIOinfo/Jappy/review/data/hiv_control/fastp_no_umi/"
path="/projects/BIOinfo/Jappy/review/data/hiv/fastp_no_umi/"
out_path_base="/projects/BIOinfo/Jappy/review/result/deduplication/nubeamdedup/"

for data_set in $datasets
do
    
    srr=${data_set%.fastq}
    
    data_set_fullpath=$(readlink -f "${path}${data_set}")
    echo ${data_set_fullpath}

    output_file_fullpath=$(readlink -f "${out_path_base}${srr}.fastq")
    ./nubeam-dedup -i ${data_set_fullpath} -o ${output_file_fullpath}

done
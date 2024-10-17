#!/bin/bash
# @Author: Pengyao Ping
# @Date:   2023-03-13 18:18:24
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-03-13 19:02:03
# ref_index=/projects/BIOinfo/Jappy/review/data/hiv_control/bowtie2_index/hiv_control_index

# datasets="processed.SRR1543964.fastq processed.SRR1543965.fastq processed.SRR1543966.fastq processed.SRR1543967.fastq processed.SRR1543968.fastq processed.SRR1543969.fastq processed.SRR1543970.fastq processed.SRR1543971.fastq"
# path="./"

# Check if at least two parameters are provided
if [ "$#" -lt 5 ]; then
    echo "Usage: $0 data_set accssion_num ref_genome ref_index out_dir"
    exit 1
fi
# Access input parameters
data_set=$1
srr=$2
ref_genome=$3
ref_index=$4
out_dir=$5

processed_data_set=${out_dir}processed.${srr}.fastq
umi_tools extract --stdin=${data_set} --bc-pattern=NNNNNNNNNNNN --log=${srr}.extract.log --stdout ${processed_data_set}

alignment_bam="${srr}.alignment.bam"

bowtie2-build ${ref_genome} $ref_index

bowtie2 -x $ref_index -U ${processed_data_set} --local -p 64 | samtools view -bS > ${alignment_bam}

sorted_alignment_bam="sorted.${srr}.alignment.bam"
samtools sort ${alignment_bam} -o ${sorted_alignment_bam}
samtools index ${sorted_alignment_bam}

dedup_bam="${srr}.dedup.bam"
dedup_log="${srr}.deduplicated"
umi_tools dedup -I ${sorted_alignment_bam} --output-stats=${dedup_log} -S ${dedup_bam}  

groups_log="${srr}.groups.tsv"
mapped_group_bam="${srr}.mapped_grouped.bam"
umi_tools group -I ${sorted_alignment_bam} --group-out=${groups_log} --output-bam -S ${mapped_group_bam}

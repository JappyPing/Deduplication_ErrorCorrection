#!/bin/bash
# @Author: Pengyao Ping
# @Date:   2023-03-13 18:18:24
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-04 14:49:40
ref_index=/projects/BIOinfo/Jappy/review/data/hiv/genome/bowtie2_index/hiv_index

path="/projects/BIOinfo/Jappy/review/data/hiv/umi_in_read/"

out_path_base="/projects/BIOinfo/Jappy/review/result/deduplication/UMI-tools/"

data_set1_fullpath=$(readlink -f "${path}SRR11207257_1.fastq")
data_set2_fullpath=$(readlink -f "${path}SRR11207257_2.fastq")

processed1=$(readlink -f "${out_path_base}processed.SRR11207257_1.fastq")
processed2=$(readlink -f "${out_path_base}processed.SRR11207257_2.fastq")
# extract
umi_tools extract --extract-method=string -I ${data_set1_fullpath} --bc-pattern=NNNNNNNNNNNN --read2-in=${data_set2_fullpath} --bc-pattern2=NNNNNNNNNNNN --stdout=${processed1} --read2-out=${processed2}

alignment_bam="SRR11207257.alignment.bam"
bowtie2 -x $ref_index -1 ${data_set1_fullpath} -2 ${data_set2_fullpath} --local -p 64 | samtools view -bS > ${alignment_bam}

sorted_alignment_bam=$(readlink -f "${out_path_base}sorted.SRR11207257.alignment.bam")
samtools sort ${alignment_bam} -o ${sorted_alignment_bam}
samtools index ${sorted_alignment_bam}

dedup_bam=$(readlink -f "${out_path_base}SRR11207257.dedup.bam")
dedup_log="SRR11207257.deduplicated"
umi_tools dedup --umi-separator=":" -I ${sorted_alignment_bam} --paired --output-stats=${dedup_log} -S ${dedup_bam}  

groups_log="SRR11207257.groups.tsv"
mapped_group_bam="SRR11207257.mapped_grouped.bam"
umi_tools group -I ${sorted_alignment_bam} --group-out=${groups_log} --output-bam -S ${mapped_group_bam}

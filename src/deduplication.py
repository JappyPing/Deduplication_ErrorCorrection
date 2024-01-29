import time
from utils import *
import subprocess
import pandas as pd
import os

def single_end_deduplication(output_dir, parallel_num):
    # all_methods = ["UMI-tools", "AmpUMI", "UMIc", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "fastp", "pRESTO", "CD-HIT-DUP", "ParDRe","Minirmd"]
    # all_methods = ["NGSReadsTreatment", "fastp", "CD-HIT-DUP", "ParDRe","Minirmd"]
    all_methods = ["BioSeqZip"]

    umi_methods = ["AmpUMI", "UMI-tools", "UMIc"] # 
    pe_methods = ["Calib", "FastUniq", "Gencore"]
    bam_methods = ["Gencore", "UMI-tools", "umitools"]
    non_umi_methods = ["CD-HIT-DUP", "FastUniq", "pRESTO", "ParDRe", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "Minirmd"]    
    datasets = ["SRR1543964", "SRR1543965", "SRR1543966","SRR1543967","SRR1543968","SRR1543969","SRR1543970","SRR1543971"]
    mis_methods = ["pRESTO", "CD-HIT-DUP", "Minirmd", "ParDRe"]
    mis_ids = [0, 1, 2, 3]

    columns = ["DataSet", "Methods", "Time", "Peak Memory"]
    time_memory_df = pd.DataFrame(columns=columns) 

    for method in all_methods:  
        # if method == "UMIc":
        #     umic_r_script = "../methods/deduplication/UMIc/R/UMIsProject2.R"
        #     command = f'Rscript {umic_r_script}'
        #     # Start the command using subprocess
        #     process = subprocess.Popen(command, shell=True)
        #     # Get the process ID
        #     pid = process.pid
        #     peak_memory_usage = 0.0
        #     try:
        #         while process.poll() is None:
        #             memory_usage = monitor_memory(pid)
        #             # Update peak memory usage if the current value is higher
        #             if memory_usage > peak_memory_usage:
        #                 peak_memory_usage = memory_usage
        #             # print(f"Elapsed Time: {elapsed_time:.2f} seconds | Memory Usage: {memory_usage:.2f} MB")
        #             time.sleep(1)
        #     finally:
        #         elapsed_time = (time.time() - start_time)/60
        #         print(f"{method}:\n")
        #         print(f"Elapsed Time: {elapsed_time:.2f} minutes | Peak Memory Usage: {peak_memory_usage:.2f} MB", flush=True)
        #         time_memory_df.loc[len(time_memory_df)] = ["all the datasets for umic", method, elapsed_time, peak_memory_usage]    
            # time_memory_df.to_csv(output_dir + "/" + method + ".time_memory_data.csv", index=False)
        # else:
        for item in datasets:
            if method in umi_methods:
                data_path = "/projects/BIOinfo/Jappy/short_read_review/data/hiv_control/fastp_umi_in_read/umi_" + item + ".fastq"  
            else:
                data_path = "/projects/BIOinfo/Jappy/short_read_review/data/hiv_control/fastp_no_umi/"+ item + ".fastq"       
            if method  not in mis_methods:
                output_path = output_dir + method + "/"
                if not os.path.exists(output_path):
                    os.makedirs(output_path)                      

            
            if method == "UMI-tools":
                ref_index="../data/hiv_control/bowtie2_index/hiv_control_index"
                method_sh = "./deduplication/UMI-tools/umi_tools_dedup.sh"
                ref_genome = "/projects/BIOinfo/steve/G38/ref/GCF_000001405.40_GRCh38.p14_genomic.fasta"
                command = f"{method_sh} {data_path} {item} {ref_genome} {ref_index} {output_path}  "      
                res = run_command(command, method, item)
                time_memory_df.loc[len(time_memory_df)] = res    
            elif method == "AmpUMI":
                ampumi = "../methods/deduplication/AmpUMI/AmpUMI.py"
                ampumi_out = output_path + item + ".fastq"
                command = f"python {ampumi} Process --fastq {data_path} --fastq_out {ampumi_out} --umi_regex 'IIIIIIIIIIII' "
                res = run_command(command, method, item)
                time_memory_df.loc[len(time_memory_df)] = res
            elif method == "NGSReadsTreatment":
                NGSReadsTre = "../methods/deduplication/NGSReadsTreatment/NgsReadsTreatment_v1.3.jar"
                command = f"java -jar {NGSReadsTre} {data_path} {parallel_num}"
                res = run_command(command, method, item)
                time_memory_df.loc[len(time_memory_df)] = res
            elif method == "Nubeam-dedup":
                Nubeam = "../methods/deduplication/nubeamdedup-master/Linux/nubeam-dedup"
                Nubeam_output = output_path + item + ".fastq" 
                command = f"{Nubeam} -i {data_path} -o {Nubeam_output}"
                res = run_command(command, method, item)
                time_memory_df.loc[len(time_memory_df)] = res
            elif method == "BioSeqZip":
                data_path = "/projects/BIOinfo/Jappy/short_read_review/data/hiv_control/fastp_no_umi/"+ item + ".fastq"
                bioseqzip = "../methods/deduplication/BioSeqZip/build/apps/bioseqzip_collapse"
                command = f"{bioseqzip} -i {data_path} -f fastq --csv-report -t {parallel_num} -m 200G -o {output_path} -b {item}"
                res = run_command(command, method, item)
                time_memory_df.loc[len(time_memory_df)] = res
            elif method == "fastp":
                data_path = "../data/hiv_control/fastp_umi_in_name/umi_" + item + ".fastq"
                fastp_output = output_path + item + ".fastq"
                command = f"fastp --dedup -A -i {data_path} -o {fastp_output} -j {item}.json -h {item}.html "
                res = run_command(command, method, item)
                time_memory_df.loc[len(time_memory_df)] = res
            # elif method == "pRESTO":
            #     for mis_n in mis_ids:
            #         output_path = output_dir + method + "/" + str(mis_n) + "/" 
            #         if not os.path.exists(output_path):
            #             os.makedirs(output_path) 
            #         presto = "/projects/BIOinfo/Jappy/short_read_review/methods/deduplication/presto-0.5.2/bin/CollapseSeq.py"
            #         log_file = item + str(mis_n) + ".log"
            #         command = f"{presto} -s {data_path} -n {mis_n} --log {log_file} --outdir {output_path}"
            elif method == "CD-HIT-DUP":
                for mis_n in mis_ids:
                    output_path = output_dir + method + "/" + str(mis_n) + "/" 
                    if not os.path.exists(output_path):
                        os.makedirs(output_path)                         
                    out_data = output_path  + item + ".fastq"
                    command = f"cd-hit-dup -i {data_path} -e {mis_n} -o {out_data}"  
                    res = run_command(command, method, item)
                    time_memory_df.loc[len(time_memory_df)] = res
            elif method == "ParDRe":
                pardre = "../methods/deduplication/ParDRe-rel2.2.5/ParDRe"
                for mis_n in mis_ids:  
                    output_path = output_dir + method + "/" + str(mis_n) + "/" 
                    if not os.path.exists(output_path):
                        os.makedirs(output_path)                                              
                    out_data = output_path  + item + ".fastq"
                    command = f"mpirun -n {parallel_num} {pardre} -i {data_path} -m {mis_n} -o {out_data}"
                    res = run_command(command, method, item)
                    time_memory_df.loc[len(time_memory_df)] = res
            elif method == "Minirmd":
                for mis_n in mis_ids:
                    output_path = output_dir + method + "/" + str(mis_n) + "/" 
                    if not os.path.exists(output_path):
                        os.makedirs(output_path)                         
                    out_data = output_path  + item + ".fastq"
                    command = f"minirmd -i {data_path} -d {mis_n} -o {out_data} -t {parallel_num}"
                    res = run_command(command, method, item)
                    time_memory_df.loc[len(time_memory_df)] = res
            #######################################################
    time_memory_df.to_csv(output_dir + "/bioseqzip.time_memory_data.csv", index=False)


def run_command(command, method, item):
    start_time = time.time()
    # Start the command using subprocess
    process = subprocess.Popen(command, shell=True)

    # Get the process ID
    pid = process.pid
    peak_memory_usage = 0.0
    try:
        while process.poll() is None:
            memory_usage = monitor_memory(pid)
            # Update peak memory usage if the current value is higher
            if memory_usage > peak_memory_usage:
                peak_memory_usage = memory_usage
            # print(f"Elapsed Time: {elapsed_time:.2f} seconds | Memory Usage: {memory_usage:.2f} MB")
            time.sleep(1)

    finally:
        elapsed_time = (time.time() - start_time)/60
        print(f"{method} on {item}: \n")
        print(f"Elapsed Time: {elapsed_time:.2f} minutes | Peak Memory Usage: {peak_memory_usage:.2f} MB", flush=True)

    return [item, method, elapsed_time, peak_memory_usage]   

    
if __name__ == "__main__": 
    output_dir = "../results/deduplication/"
    parallel_num = 64
    single_end_deduplication(output_dir, parallel_num)
import time
from utils import *
import subprocess
import pandas as pd
import os

def single_end_deduplication(output_dir, parallel_num):
    # all_methods = ["UMI-tools", "AmpUMI", "UMIc", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "fastp", "pRESTO", "CD-HIT-DUP", "ParDRe","Minirmd"]
    all_methods = ["pRESTO"]

    # umi_methods = ["AmpUMI", "BioSeqZip", "fastp", "UMI-tools", "UMIc"] # 
    # pe_methods = ["Calib", "FastUniq", "Gencore"]
    # bam_methods = ["Gencore", "UMI-tools", "umitools"]
    # non_umi_methods = ["CD-HIT-DUP", "FastUniq", "pRESTO", "ParDRe", "NGSReadsTreatment", "Nubeam-dedup", "BioSeqZip", "Minirmd"]    
    datasets = ["SRR1543964", "SRR1543965", "SRR1543966","SRR1543967","SRR1543968","SRR1543969","SRR1543970","SRR1543971"]
    # mis_methods = ["pRESTO", "CD-HIT-DUP", "Minirmd", "ParDRe"]
    mis_ids = [0, 1, 2, 3]

    columns = ["DataSet", "Methods", "Time", "Peak Memory"]
    time_memory_df = pd.DataFrame(columns=columns) 

    for item in datasets:
        data_path = "/projects/BIOinfo/Jappy/short_read_review/data/hiv_control/fastp_no_umi/"+ item + ".fastq"       

        start_time = time.time()

        for mis_n in mis_ids:
            output_path ="/projects/BIOinfo/Jappy/short_read_review/results/deduplication/pRESTO/" + str(mis_n) 
            # if not os.path.exists(output_path):
            #     os.makedirs(output_path) 
            log_file = item + str(mis_n) + ".log"
            command = f"CollapseSeq.py -s {data_path} -n {mis_n} --log {log_file} --outdir {output_path}"
        #######################################################
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
                print(f"pRESTO on {item}: \n")
                print(f"Elapsed Time: {elapsed_time:.2f} minutes | Peak Memory Usage: {peak_memory_usage:.2f} MB", flush=True)
                time_memory_df.loc[len(time_memory_df)] = [item, "pRESTO", elapsed_time, peak_memory_usage]   
    time_memory_df.to_csv(output_dir + "/time_memory_data.csv", index=False)
if __name__ == "__main__": 
    output_dir = "/projects/BIOinfo/Jappy/short_read_review/results/deduplication/"
    parallel_num = 64
    single_end_deduplication(output_dir, parallel_num)
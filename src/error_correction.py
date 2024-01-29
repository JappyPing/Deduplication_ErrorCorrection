import os
import time
# from datetime import datetime
from utils import *
import subprocess
import pandas as pd

def correction(method, raw_dir, correct_base_dir, parallel_n):
    columns = ["DataSet", "Methods", "Time", "Peak Memory"]
    time_memory_df = pd.DataFrame(columns=columns) 

    file_names = os.listdir(raw_dir)
    correct_dir = correct_base_dir + method + "/"
    if not os.path.exists(correct_dir):
        os.makedirs(correct_dir)

    for file_name in file_names:
        start_time = time.time()
        raw = os.path.join(raw_dir, file_name)
        correct = os.path.join(correct_dir, file_name)
        genome_len = 10000
        if method == "BFC":
            bin_method = "../methods/errorCorrection/bfc"
            k = 30
            # os.system("%s -t26 -k %s %s > %s" % (bin_method, k, raw, correct))
            command = f"{bin_method} -t {parallel_n} -k {k} {raw} > {correct}"
        elif method == "Coral":
            bin_method = "../methods/errorCorrection/coral"
            if "SRR1543965" in file_name:
                k = 28
            else:
                k = 30
            # os.system("%s -fq %s -p 20 -o %s -k %s" % (bin_method, raw, correct, k))
            command = f"{bin_method} -fq {raw} -p {parallel_n} -o {correct} -k {k}"
        elif method == "Fiona":
            bin_method = "../methods/errorCorrection/fiona"
            # os.system("%s -nt 20 -g %s %s %s" % (bin_method, genome_len, raw, correct))
            command = f"{bin_method} -nt {parallel_n} -g {genome_len} {raw} {correct}"
        elif method == "Lighter":
            bin_method = "../methods/errorCorrection/lighter"
            if "SRR1543967" in file_name or "SRR1543970" in file_name:
                k = 28
            elif "SRR1543968" in file_name:
                k = 26
            else:
                k = 30
            # os.system("%s -t 20 -r %s -K %s %s -od %s" % (bin_method, raw, k, genome_len, correct_dir)) 
            command = f"{bin_method} -t {parallel_n} -r {raw} -K {k} {genome_len} -od {correct_dir}"  
        # elif method == "musket":
        #     k = 30
        #     os.system("%s -k %s 134217728 -o %s -p 20 %s" % (bin_method, k, correct,raw))        
        elif method == "Pollux":
            bin_method = "../methods/errorCorrection/pollux"
            k = 20
            # os.system("%s -i %s -o %s -k %s" % (bin_method, raw, correct_dir, k))   
            command = f"{bin_method} -i {raw} -o {correct_dir} -k {k}"
            # Start the command using subprocess
            process = subprocess.Popen(command, shell=True)

            # Get the process ID
            pid = process.pid
            peak_memory_usage = 0.0
            try:
                while process.poll() is None:
                    elapsed_time = time.time() - start_time
                    memory_usage = monitor_memory(pid)
                    
                    # Update peak memory usage if the current value is higher
                    if memory_usage > peak_memory_usage:
                        peak_memory_usage = memory_usage
                    
                    # print(f"Elapsed Time: {elapsed_time:.2f} seconds | Memory Usage: {memory_usage:.2f} MB")
                    
                    time.sleep(1)

            finally:
                elapsed_time = (time.time() - start_time)/60
                print(f"Elapsed Time: {elapsed_time:.2f} minutes | Peak Memory Usage: {peak_memory_usage:.2f} MB", flush=True)
                time_memory_df.loc[len(time_memory_df)] = [file_name, method, elapsed_time, peak_memory_usage]

            os.system("cat %s %s > %s" % (correct_dir + file_name + ".corrected", correct_dir + file_name + ".low", correct_dir + file_name))
        elif method == "RACER":
            bin_method = "../methods/errorCorrection/RACER"
            # os.system("%s %s %s %s" % (bin_method, raw, correct, genome_len))   
            command = f"{bin_method} {raw} {correct} {genome_len}"
        elif method == "Karect":
            bin_method = "../methods/errorCorrection/karect"
            command = f"{bin_method} -correct -threads={parallel_n} -matchtype=edit -celltype=haploid  -minoverlap=120 -kmer=14 -kmererrors=2 -inputfile={raw} -resultdir={correct_dir}"
        #     os.system("%s -correct -threads=26 -matchtype=edit -celltype=haploid  -minoverlap=120 -kmer=14 -kmererrors=2 -inputfile=%s -resultdir=%s" % (bin_method, raw, correct_dir)) 
        elif method == "Care":
            bin_method = "../methods/errorCorrection/care-cpu"
            coverage = 60
            # os.system("%s -p -i %s -c %s -o %s -d %s -m 27G -t 20 --pairmode SE -h 48" % (bin_method, raw, coverage, file_name, correct_dir)) 
            command = f"{bin_method} -p -i {raw} -c {coverage} -o {file_name} -d {correct_dir} -m 27G -t {parallel_n} --pairmode SE -h 48"
        elif method == "Bcool":
            base_name = file_name.split(".fastq")[0]
            raw_fa_dir = correct_dir + "raw_fa/" 
            if not os.path.exists(raw_fa_dir):
                os.makedirs(raw_fa_dir)           
            raw_fa = raw_fa_dir + base_name + ".fasta" 
            # covert_fq2fa(raw, raw_fa)   

            bcool_correct_dir = os.path.join(correct_dir, file_name + "/input/")
            if not os.path.exists(bcool_correct_dir):
                os.makedirs(bcool_correct_dir)
            # print(raw_fa)

            start_time = time.time()
            # os.system("bcool -t 60 -u %s -o %s" % (raw_fa, bcool_correct_dir))
            command = f"bcool -t {parallel_n} -u {raw_fa} -o {bcool_correct_dir}"

        if method != "Pollux":
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
                print(f"Elapsed Time: {elapsed_time:.2f} minutes | Peak Memory Usage: {peak_memory_usage:.2f} MB", flush=True)
                time_memory_df.loc[len(time_memory_df)] = [file_name, method, elapsed_time, peak_memory_usage]

    return time_memory_df

def main(raw_dir, correct_base_dir):
    columns = ["DataSet", "Methods", "Time", "Peak Memory"]
    time_memory_df = pd.DataFrame(columns=columns) 

    methods = ["Bcool", "Care", "Coral",  "Fiona", "Lighter", "RACER", "BFC", "Karect"]
    # methods = ["Pollux"]
    parallel_n = 64

    for method in methods:
        print("#############################################################")
        print(method)
        tm_result = correction(method, raw_dir, correct_base_dir, parallel_n)

        time_memory_df = pd.concat([time_memory_df, tm_result])
        # time_memory_df.loc[len(time_memory_df)] = tm_result
        tm_result.to_csv(correct_base_dir + method + "/time_memory_data.csv", index=False)
        print("#############################################################")

    # Save time_memory_df to a file (e.g., CSV)
    time_memory_df.to_csv(correct_base_dir + "/time_memory_data.csv", index=False)

if __name__ == "__main__":
    hiv_control_raw_dir = "../data/hiv_control/fastp_no_umi/"
    correct_base_dir = "../results/error_correction/"
    main(hiv_control_raw_dir, correct_base_dir)

    hiv_raw_dir = "../data/hiv/fastp_no_umi/"





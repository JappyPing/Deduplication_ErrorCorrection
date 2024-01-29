# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-05-10 10:01:27
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-11 16:53:11
import os

def main():
    mis_matches_lst = ["1", "2", "3"]
    raw_dir = "/projects/BIOinfo/Jappy/review/data/hiv_control/fastp_no_umi/"
    file_names = os.listdir(raw_dir)
    out_dir = "/projects/BIOinfo/Jappy/review/result/deduplication/pRESTO/"
    for mis in mis_matches_lst:
        cur_out_dir = out_dir + mis + "/"
        if not os.path.exists(cur_out_dir):
            os.makedirs(cur_out_dir)

        for f_n in file_names:
            f_n_path = raw_dir + f_n
            # print(f_n_path)
            log_file = f_n + mis + ".log"
            presto_commands = f"CollapseSeq.py -s {f_n_path} -n {int(mis)} --log {log_file} --outdir {cur_out_dir}"
            try:
                os.system(presto_commands)
            except OSError as e:
                print(f"Error executing command: {e}")


if __name__ == "__main__":
    main()

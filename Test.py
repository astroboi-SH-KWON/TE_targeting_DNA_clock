import time
import os
from Bio import SeqIO
import multiprocessing as mp
import numpy as np
import platform

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = os.getcwd() + "/"
SYSTEM_NM = platform.system()
if SYSTEM_NM == 'Linux':
    # REAL
    REF_DIR = "../hg38/"
    DFAM_ANNO = "./input/hg38_dfam.hits"
else:
    # DEV
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"
    DFAM_ANNO = "D:/000_WORK/ParkJiHye/20200914/hg38_dfam.hits"

PROJECT_NAME = WORK_DIR.split("/")[-2]
FILTERED_CDS_INFO = "filtered_hg38_refFlat.txt"
multi_processing_1_FILE = "ClinVar_hg38_result.txt"

# name, pam_seq, len_spacer, win_size_arr
INIT = [
    ['TE_trgt', 'NGG', 20, [10, 10]]
]

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
############### end setting env #################


def test():
    try:
        os.remove(DFAM_ANNO + "_top_500.txt")
        os.remove(DFAM_ANNO + "_tail_500.txt")
    except Exception as err:
        print(str(err))

    with open(DFAM_ANNO, 'r') as in_f:
        with open(DFAM_ANNO + "_top_500.txt", 'a') as out_top_500_f:
            with open(DFAM_ANNO + "_tail_500.txt", 'a') as out_tail_500_f:
                frst_line = in_f.readline()
                out_top_500_f.write(frst_line)
                out_tail_500_f.write(frst_line)
                chr_nm = ""
                idx = 0
                tmp_list = []
                while True:
                    tmp_line = in_f.readline()
                    if tmp_line == '':
                        break
                    if chr_nm != tmp_line.split('\t')[0]:
                        print(chr_nm, ": chr_nm")
                        chr_nm = tmp_line.split('\t')[0]
                        idx = 0
                        for tp_line in tmp_list:
                            out_tail_500_f.write(tp_line)

                    tmp_list.append(tmp_line)
                    if len(tmp_list) > 500:
                        tmp_list = tmp_list[1:]
                    if idx < 500:
                        out_top_500_f.write(tmp_line)
                    idx += 1









if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    test()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
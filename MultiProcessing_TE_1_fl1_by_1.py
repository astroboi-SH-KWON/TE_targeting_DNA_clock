import time
import os
from Bio import SeqIO
import multiprocessing as mp
import numpy as np
import platform
import math

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = os.getcwd() + "/"
SYSTEM_NM = platform.system()
if SYSTEM_NM == 'Linux':
    # REAL
    REF_DIR = "../hg38/"
    DFAM_ANNO = "./input/hg38_dfam.nrph.hits"
else:
    # DEV
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"
    DFAM_ANNO = "D:/000_WORK/ParkJiHye/20200914/hg38_dfam.nrph.hits"  # 55

PROJECT_NAME = WORK_DIR.split("/")[-2]
############### end setting env #################

def main_TE_1_fl1_by_1():
    util = Util.Utils()
    header = ['sequence', '#duple', '#trnscprt', 'chromosome:23bp(spacer + PAM)index range:strand:transcription:fam_name']

    fl_sources = util.get_files_from_dir(WORK_DIR + "output/TE_trgt_*.txt")
    for trgt_fl in fl_sources:
        dfam_list = util.read_csv_ignore_N_line(trgt_fl, '\t')
        result_dict = {}
        for dfam_arr in dfam_list:
            chr_nm = dfam_arr[0]
            tot_seq = dfam_arr[1]
            fam_nm = dfam_arr[2]
            index_rng = dfam_arr[3]
            strand = dfam_arr[4]
            trns_flag = dfam_arr[5]

            res_key = tot_seq
            if res_key in result_dict:
                result_dict[res_key].add(chr_nm + ":" + index_rng + ":" + strand + ":" + trns_flag + ":" + fam_nm)
            else:
                tmp_set = set()
                tmp_set.add(chr_nm + ":" + index_rng + ":" + strand + ":" + trns_flag + ":" + fam_nm)
                result_dict.update({res_key: tmp_set})

        result_list = []
        for res_key, val_set in result_dict.items():
            tmp_str = ""
            cnt_trpt = 0
            for tmp_val in val_set:
                if 'True' in tmp_val:
                    cnt_trpt += 1
                tmp_str += tmp_val + ", "
            result_list.append([res_key, len(val_set), cnt_trpt, tmp_str])

        util.make_csv(
            trgt_fl[:trgt_fl.index("TE_trgt_")] + "loop/" + trgt_fl[trgt_fl.index("TE_trgt_"):].replace(".txt", "") + ".fl1_by_1", header, result_list, 0, '\t')

def main_TE_1_fl_n_by_1():
    n = 10
    util = Util.Utils()
    header = ['sequence', '#duple', '#trnscprt',
              'chromosome:23bp(spacer + PAM)index range:strand:transcription:fam_name']

    fl_sources = util.get_files_from_dir(WORK_DIR + "output/loop/TE_trgt_*.fl1_by_1")

    last_epoch = len(fl_sources) % n
    epoch = math.ceil(len(fl_sources)/n)
    trgt_fl = fl_sources[0]
    for i in range(epoch):
        dfam_list = []
        if (i + 1) == epoch:
            if last_epoch == 0:
                for j in range(i * n, (i + 1) * n):
                    dfam_list.extend(util.read_csv_ignore_N_line(fl_sources[j], '\t'))
            else:
                for j in range(i * n, i * n + last_epoch):
                    dfam_list.extend(util.read_csv_ignore_N_line(fl_sources[j], '\t'))
        else:
            for j in range(i*n, (i+1)*n):
                dfam_list.extend(util.read_csv_ignore_N_line(fl_sources[j], '\t'))

        result_dict = {}
        for dfam_arr in dfam_list:
            tot_seq = dfam_arr[0]
            if 'N' in tot_seq:
                continue

            res_key = tot_seq
            if res_key in result_dict:
                result_dict[res_key].update(dfam_arr[3].replace(" ", "").split(',')[:-1])
            else:
                tmp_set = set(dfam_arr[3].replace(" ", "").split(',')[:-1])
                result_dict.update({res_key: tmp_set})

        result_list = []
        for res_key, val_set in result_dict.items():
            tmp_str = ""
            cnt_trpt = 0
            for tmp_val in val_set:
                if 'True' in tmp_val:
                    cnt_trpt += 1
                tmp_str += tmp_val + ", "
            result_list.append([res_key, len(val_set), cnt_trpt, tmp_str])

        util.make_csv(trgt_fl[:trgt_fl.index("TE_trgt_")] + "TE_trgt_cyc1_" + str(i) + "_fl" + str(n) + "_by_1", header,
                      result_list, 0, '\t')

def main_TE_1_fl_n_by_1_w_array():
    num_arr = [[0,1,2],[3,4,5],[6,7,8],[9,10],[11,12]]
    util = Util.Utils()
    header = ['sequence', '#duple', '#trnscprt',
              'chromosome:23bp(spacer + PAM)index range:strand:transcription:fam_name']

    fl_nm_f = WORK_DIR + "output/loop/TE_trgt_cyc"
    fl_nm_b = "_fl2_by_1"
    for i in range(len(num_arr)):
        dfam_list = []
        for j in num_arr[i]:
            dfam_list.extend(util.read_csv_ignore_N_line(fl_nm_f + "2_" + str(j) + fl_nm_b, '\t'))

        result_dict = {}
        for dfam_arr in dfam_list:
            tot_seq = dfam_arr[0]
            if 'N' in tot_seq:
                continue

            res_key = tot_seq
            if res_key in result_dict:
                result_dict[res_key].update(dfam_arr[3].replace(" ", "").split(',')[:-1])
            else:
                tmp_set = set(dfam_arr[3].replace(" ", "").split(',')[:-1])
                result_dict.update({res_key: tmp_set})

        result_list = []
        for res_key, val_set in result_dict.items():
            tmp_str = ""
            cnt_trpt = 0
            for tmp_val in val_set:
                if 'True' in tmp_val:
                    cnt_trpt += 1
                tmp_str += tmp_val + ", "
            result_list.append([res_key, len(val_set), cnt_trpt, tmp_str])

        util.make_csv(fl_nm_f + "3_" + str(i) + "_fln" + "_by_1", header, result_list, 0, '\t')

def main_TE_1_fl_n_by_1_right_away():
    util = Util.Utils()
    header = ['sequence', '#duple', '#trnscprt',
              'chromosome:23bp(spacer + PAM)index range:strand:transcription:fam_name']

    fl_nm_f = WORK_DIR + "output/loop/TE_trgt_cyc"
    fl_nm_b = "_fln_by_1"

    result_dict = {}
    tm_arr = [[0,3,4]]
    for i_a in range(len(tm_arr)):
        for i in tm_arr[i_a]:
            print(fl_nm_f + "3_" + str(i) + fl_nm_b)
            with open(fl_nm_f + "3_" + str(i) + fl_nm_b) as f:
                print(f.readline())
                while True:
                    tmp_line = f.readline().replace("\n", "")
                    if tmp_line == '':
                        break
                    dfam_arr = tmp_line.split('\t')
                    tot_seq = dfam_arr[0]
                    res_key = tot_seq
                    if res_key in result_dict:
                        result_dict[res_key].update(dfam_arr[3].replace(" ", "").split(',')[:-1])
                    else:
                        tmp_set = set(dfam_arr[3].replace(" ", "").split(',')[:-1])
                        result_dict.update({res_key: tmp_set})

        result_list = []
        for res_key, val_set in result_dict.items():
            tmp_str = ""
            cnt_trpt = 0
            for tmp_val in val_set:
                if 'True' in tmp_val:
                    cnt_trpt += 1
                tmp_str += tmp_val + ", "
            result_list.append([res_key, len(val_set), cnt_trpt, tmp_str])

        util.make_csv(fl_nm_f + "4_" + str(i_a) + "_fln_by_1", header, result_list, 0, '\t')

def split_TE_1_fl_n_by_1_right_away():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    header = ['sequence', '#duple', '#trnscprt',
              'chromosome:23bp(spacer + PAM)index range:strand:transcription:fam_name']

    fl_nm_f = WORK_DIR + "output/loop/TE_trgt_cyc"
    cyc_num = 3
    fl_nm_b = "_fln_by_1"

    res_f_num = 0

    # tm_arr = [[1, 2, 3, 5], [0, 4, 6, 7]]
    # tm_arr = [[0, 3], [1, 2]]
    tm_arr = [[0, 2]]
    for i_a in range(len(tm_arr)):
        result_dict = {}
        for i in tm_arr[i_a]:
            print(fl_nm_f + str(cyc_num) + "_" + str(i) + fl_nm_b)
            with open(fl_nm_f + str(cyc_num) + "_" + str(i) + fl_nm_b) as f:
                print(f.readline())
                while True:
                    tmp_line = f.readline().replace("\n", "")
                    if tmp_line == '':
                        break
                    dfam_arr = tmp_line.split('\t')
                    tot_seq = dfam_arr[0]
                    res_key = tot_seq
                    if res_key in result_dict:
                        result_dict[res_key].update(dfam_arr[3].replace(" ", "").split(',')[:-1])
                    else:
                        tmp_set = set(dfam_arr[3].replace(" ", "").split(',')[:-1])
                        result_dict.update({res_key: tmp_set})

        result0_list = []
        result1_list = []
        for res_key, val_set in result_dict.items():
            tmp_str = ""
            cnt_trpt = 0
            for tmp_val in val_set:
                if 'True' in tmp_val:
                    cnt_trpt += 1
                tmp_str += tmp_val + ", "
            if len(val_set) > 1:
                result0_list.append([res_key, len(val_set), cnt_trpt, tmp_str])
            else:
                result1_list.append([res_key, len(val_set), cnt_trpt, tmp_str])

        result_dict.clear()
        sorted_result0_list = logic_prep.sort_list_by_ele(result0_list, 1)
        result0_list.clear()

        util.make_csv(fl_nm_f + str(cyc_num + 1) + "_" + str(res_f_num) + "_fln_by_1", header, sorted_result0_list, 0, '\t')
        res_f_num += 1
        util.make_csv(fl_nm_f + str(cyc_num + 1) + "_" + str(res_f_num) + "_fln_by_1", header, result1_list, 0, '\t')
        res_f_num += 1
        sorted_result0_list.clear()
        result1_list.clear()

def multi_processing_TE_1_fl_n_by_1():
    num_proc = 7

    util = Util.Utils()
    header = ['sequence', '#duple', '#trnscprt',
              'chromosome:23bp(spacer + PAM)index range:strand:transcription:fam_name']

    # fl_sources = util.get_files_from_dir(WORK_DIR + "output/loop/TE_trgt_*.fl1_by_1")

    num_arr = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10], [11, 12]]

    for i in range(len(num_arr)):
        dfam_list = []
        for j in num_arr[i]:
            print(j)

def split_file():
    big_f = WORK_DIR + "output/loop/TE_trgt_cyc4_0_fln_by_1"
    num_split = 30
    max_row = 600

    with open(big_f) as input_f:
        header = input_f.readline()
        for num in range(num_split):
            with open(WORK_DIR + "output/loop/TE_trgt_" + str(num), 'w') as out_f:
                out_f.write(header)
                cnt = 0
                for tmp_line in input_f:
                    cnt += 1
                    out_f.write(tmp_line)
                    if cnt == max_row:
                        break
            if num == 0:
                max_row = 2000
            elif num == 1:
                max_row = 10000
            elif num == 2:
                max_row = 30000
            elif num == 3:
                max_row = 80000
            elif num == 4:
                max_row = 170000
            elif num == 5:
                max_row = 300000
            elif num == 6:
                max_row = 450000
            elif num < 10:
                max_row = 600000
            else:
                max_row = 800000


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # main_TE_1_fl1_by_1()
    # main_TE_1_fl_n_by_1()
    # main_TE_1_fl_n_by_1_w_array()
    # main_TE_1_fl_n_by_1_right_away()
    # split_TE_1_fl_n_by_1_right_away()
    split_file()
    # multi_processing_TE_1_fl_n_by_1()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))

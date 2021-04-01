import time
import os
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
else:
    # DEV
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"
    WORK_DIR = "D:/000_WORK/ParkJiHye/20210330/WORK_DIR/"

PROJECT_NAME = WORK_DIR.split("/")[-2]
FILTERED_CDS_INFO = "filtered_hg38_refFlat.txt"
TE_info_fl = "Human Full Genome_TandemRepeat_TRD_20210330.txt"

IN = "input/"
OU = "output/"
os.makedirs(WORK_DIR + IN, exist_ok=True)
os.makedirs(WORK_DIR + OU, exist_ok=True)

# name, pam_seq, len_spacer, win_size_arr
INIT = ['TE_trgt', 'NGG', 20, [10, 10]]

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)

CDS_INF = Util.Utils().read_csv_ignore_N_line(WORK_DIR + IN + FILTERED_CDS_INFO, "\t")
############### end setting env #################


"""
    + : 5'- spacer(20nt) - PAM -3'
    - : 3'- PAM - spacer(20nt) -5'
"""
def get_guide_seq_idx_strnd_trns_flg(te_inf_arr, cds_info, result_list):
    logic = Logic.Logics()

    pam_rule = INIT[1]
    len_pam = len(pam_rule)
    len_spacer = INIT[2]
    len_f_win = INIT[3][0]
    len_b_win = INIT[3][1]

    chr_nm = te_inf_arr[15]
    regn_seq_lft = te_inf_arr[16].upper()
    trgt_p_seq = te_inf_arr[18].upper()
    trgt_m_seq = ""
    try:
        trgt_m_seq += logic.make_complement_string(trgt_p_seq)
    except Exception as err:
        print("trgt_m_seq", trgt_m_seq)
        print(te_inf_arr, err)
    regn_seq_rgt = te_inf_arr[19].upper()
    fam_nm = str(te_inf_arr[0] + "_" + te_inf_arr[20])
    env_st = int(te_inf_arr[1])
    # last_idx = int(te_inf_arr[2])

    p_trgt_seq_f = regn_seq_lft[- len_f_win:]
    p_trgt_seq_b = regn_seq_rgt[: len_b_win]

    m_trgt_seq_f = ""
    try:
        m_trgt_seq_f += logic.make_complement_string(regn_seq_rgt)[: len_f_win]
    except Exception as e:
        print("m_trgt_seq_f", m_trgt_seq_f)
        print(te_inf_arr, e)

    m_trgt_seq_b = ""
    try:
        m_trgt_seq_b += logic.make_complement_string(regn_seq_lft)[- len_b_win:]
    except Exception as e:
        print("m_trgt_seq_b", m_trgt_seq_b)
        print(te_inf_arr, e)

    for i in range(len(trgt_p_seq) - len_spacer - len_pam + 1):
        p_pam = trgt_p_seq[i + len_spacer: i + len_spacer + len_pam]
        m_pam = trgt_m_seq[i: i + len_pam]

        real_pos_st = env_st + i
        real_pos_en = env_st + i + len_spacer + len_pam

        # check + strand
        if logic.match(0, p_pam, pam_rule):
            spacer = trgt_p_seq[i: i + len_spacer]

            tmp_idx = i - len_b_win
            if tmp_idx < 0:
                tmp_idx = 0
            f_win_seq = trgt_p_seq[tmp_idx: i]
            if len(f_win_seq) < len_f_win:
                f_win_seq = p_trgt_seq_f[-(len_f_win - len(f_win_seq)):] + f_win_seq

            b_win_seq = trgt_p_seq[i + len_spacer + len_pam: i + len_spacer + len_pam + len_b_win]
            if len(b_win_seq) < len_b_win:
                b_win_seq += p_trgt_seq_b[:len_b_win - len(b_win_seq)]

            tot_seq = f_win_seq + spacer + p_pam + b_win_seq

            if not logic.is_N_in_seq(tot_seq):
                trns_flag = False
                for cds_arr in cds_info:
                    if chr_nm != cds_arr[2]:
                        continue

                    # gen_sym = cds_arr[0]
                    # nm_id = cds_arr[1]
                    trns_st = int(cds_arr[4])
                    trns_en = int(cds_arr[5])
                    if trns_st < real_pos_st and real_pos_en < trns_en:
                        trns_flag = True
                        break

                result_list.append([chr_nm, tot_seq, fam_nm, str(real_pos_st) + "-" + str(real_pos_en), '+', trns_flag])

        # check - strand
        if logic.match(0, m_pam, pam_rule[::-1]):
            spacer = trgt_m_seq[i + len_pam: i + len_pam + len_spacer]

            f_win_seq = trgt_m_seq[i + len_pam + len_spacer: i + len_pam + len_spacer + len_f_win]
            if len(f_win_seq) < len_f_win:
                f_win_seq += m_trgt_seq_f[: len_f_win - len(f_win_seq)]

            tmp_idx = i - len_b_win
            if tmp_idx < 0:
                tmp_idx = 0
            b_win_seq = trgt_m_seq[tmp_idx: i]
            if len(b_win_seq) < len_b_win:
                b_win_seq = m_trgt_seq_b[- (len_b_win - len(b_win_seq)):] + b_win_seq

            tot_seq = (b_win_seq + m_pam + spacer + f_win_seq)[::-1]


            if not logic.is_N_in_seq(tot_seq):
                trns_flag = False
                for cds_arr in cds_info:
                    if chr_nm != cds_arr[2]:
                        continue

                    # gen_sym = cds_arr[0]
                    # nm_id = cds_arr[1]
                    trns_st = int(cds_arr[4])
                    trns_en = int(cds_arr[5])
                    if trns_st < real_pos_st and real_pos_en < trns_en:
                        trns_flag = True
                        break

                result_list.append([chr_nm, tot_seq, fam_nm, str(real_pos_st) + "-" + str(real_pos_en), '-', trns_flag])


def start_multi_processing(te_info_list):
    result_list = []
    print("st_multi_processing")

    for te_inf_arr in te_info_list:
        get_guide_seq_idx_strnd_trns_flg(te_inf_arr, CDS_INF, result_list)
    print("DONE_multi_processing")
    return result_list


def main_by_list():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    sources = util.get_files_from_dir(WORK_DIR + IN + TE_info_fl.replace(".txt", "") + "/Genome_TandemRepeat_TRD*.txt")

    header = ['chr', 'tot_seq', 'fam_nm', 'index', 'strand', 'trns_flag']
    for i in range(len(sources)):
        te_inf_list = util.read_csv_ignore_N_line(sources[i], "\t", 0)

        splited_te_inf_list = np.array_split(te_inf_list, MULTI_CNT)

        print("platform.system() : ", SYSTEM_NM)
        print("total cpu_count : ", str(TOTAL_CPU))
        print("will use : ", str(MULTI_CNT))
        pool = mp.Pool(processes=MULTI_CNT)

        pool_list = pool.map(start_multi_processing, splited_te_inf_list)
        pool.close()
        splited_te_inf_list[:] = []

        result_list = logic_prep.merge_multi_list(pool_list)
        pool_list.clear()

        util.make_csv(WORK_DIR + "output/TE_trgt_20210330_" + str(i) + ".txt", header, result_list, 0, '\t')
        result_list.clear()


def main_by_list_w_filenames():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    # file_num_list = []
    # for j in range(964):
    #     file_num_list.append(j)

    header = ['chr', 'tot_seq', 'fam_nm', 'index', 'strand', 'trns_flag']
    # for i in file_num_list:
    for i in range(964):
        path = WORK_DIR + IN + TE_info_fl.replace(".txt", "") + "/Genome_TandemRepeat_TRD_" + str(i) + ".txt"
        te_inf_list = util.read_csv_ignore_N_line(path, "\t", 0)

        splited_te_inf_list = np.array_split(te_inf_list, MULTI_CNT)

        print("platform.system() : ", SYSTEM_NM)
        print("total cpu_count : ", str(TOTAL_CPU))
        print("will use : ", str(MULTI_CNT))
        pool = mp.Pool(processes=MULTI_CNT)

        pool_list = pool.map(start_multi_processing, splited_te_inf_list)
        pool.close()
        splited_te_inf_list[:] = []

        result_list = logic_prep.merge_multi_list(pool_list)
        pool_list.clear()

        util.make_csv(WORK_DIR + "output2/TE_trgt_20210330_" + str(i) + ".txt", header, result_list, 0, '\t')
        result_list.clear()


def split_big_file(num_row):
    print("split_big_file", num_row)
    util = Util.Utils()

    init_split_file = {'big_file_path': WORK_DIR + IN + TE_info_fl.replace(".txt", "") + "_wout_header.txt"
                        , 'num_row': num_row
                        , 'splited_files_dir': WORK_DIR + IN + TE_info_fl.replace(".txt", "") + "/"
                        , 'output_file_nm': "Genome_TandemRepeat_TRD"
                        , 'output_file_ext': '.txt'
                       }

    util.split_big_file_by_row(init_split_file)


def test():
    # with open(WORK_DIR + IN + TE_info_fl, 'r') as f:
    #     print(f.readline())
    with open(WORK_DIR + IN + TE_info_fl.replace(".txt", "") + "_wout_header.txt", 'r') as f:

        while True:
            tmp_line = f.readline().replace("\n", "")
            if tmp_line == "":
                break

            te_inf_arr = tmp_line.split("\t")
            # tot_seq = str(te_inf_arr).upper()
            regn_seq_lft = te_inf_arr[16].upper()
            trgt_p_seq = te_inf_arr[18].upper()
            regn_seq_rgt = te_inf_arr[19].upper()
            tot_seq = regn_seq_lft + trgt_p_seq + regn_seq_rgt

            if "R" in tot_seq:
                print(tmp_line)


def remove_head(ignr_N_line):
    with open(WORK_DIR + IN + TE_info_fl, 'r') as f_in:
        for i in range(ignr_N_line):
            f_in.readline()

        with open(WORK_DIR + IN + TE_info_fl.replace(".txt", "") + "_wout_header.txt", "w") as f_ou:
            while True:
                tmp_line = f_in.readline().replace("\n", "")
                if tmp_line == "":
                    break

                f_ou.write(tmp_line + "\n")


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # remove_head(1)
    # split_big_file(1000)
    # main_by_list()
    main_by_list_w_filenames()
    # test()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))

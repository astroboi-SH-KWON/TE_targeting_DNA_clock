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
    trgt_m_seq = logic.make_complement_string(trgt_p_seq)
    regn_seq_rgt = te_inf_arr[19].upper()
    fam_nm = str(te_inf_arr[0] + "_" + te_inf_arr[20])
    env_st = int(te_inf_arr[1])
    last_idx = int(te_inf_arr[2])

    p_trgt_seq_f = regn_seq_lft[- len_f_win:]
    p_trgt_seq_b = regn_seq_rgt[: len_b_win]

    m_trgt_seq_f = logic.make_complement_string(regn_seq_rgt)[: len_f_win]
    m_trgt_seq_b = logic.make_complement_string(regn_seq_lft)[ - len_b_win:]

    for i in range(len(trgt_p_seq) - len_spacer - len_pam + 1):
        p_pam = trgt_p_seq[i + len_spacer: i + len_spacer + len_pam]
        m_pam = trgt_m_seq[i: i + len_pam]

        real_pos_st = env_st + i
        real_pos_en = env_st + i + len_spacer + len_pam

        # check + strand
        if logic.match(0, p_pam, pam_rule):
            spacer = trgt_p_seq[i: i + len_spacer]

            f_win_seq = trgt_p_seq[i - len_f_win: i]
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

                    gen_sym = cds_arr[0]
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

            b_win_seq = trgt_m_seq[i - len_b_win: i]
            if len(b_win_seq) < len_b_win:
                b_win_seq = m_trgt_seq_b[- (len_b_win - len(b_win_seq)):] + b_win_seq

            tot_seq = (b_win_seq + m_pam + spacer + f_win_seq)[::-1]

            if not logic.is_N_in_seq(tot_seq):
                trns_flag = False
                for cds_arr in cds_info:
                    if chr_nm != cds_arr[2]:
                        continue

                    gen_sym = cds_arr[0]
                    # nm_id = cds_arr[1]
                    trns_st = int(cds_arr[4])
                    trns_en = int(cds_arr[5])
                    if trns_st < real_pos_st and real_pos_en < trns_en:
                        trns_flag = True
                        break

                result_list.append([chr_nm, tot_seq, fam_nm, str(real_pos_st) + "-" + str(real_pos_en), '-', trns_flag])


def main():
    util = Util.Utils()
    with open(WORK_DIR + IN + TE_info_fl, 'r') as f:
        cds_inf = util.read_csv_ignore_N_line(WORK_DIR + IN + FILTERED_CDS_INFO, "\t")
        result_list = []
        print(f.readline().replace("\n", ""))
        cnt = 0
        while True:
            cnt += 1
            tmp_line = f.readline().replace("\n", "")
            if tmp_line == "":
                break
            if SYSTEM_NM != 'Linux':
                # DEV
                if cnt == 1000:
                    break
            te_info_arr = tmp_line.split("\t")

            get_guide_seq_idx_strnd_trns_flg(te_info_arr, cds_inf, result_list)

        header = ['chr', 'tot_seq', 'fam_nm', 'index', 'strand', 'trns_flag']
        util.make_csv(WORK_DIR + OU + "TE_trgt_20210330.txt", header, result_list, 0, '\t')


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    main()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
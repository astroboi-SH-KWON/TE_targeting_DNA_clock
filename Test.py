import time
import os
# from Bio import SeqIO
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
    DFAM_ANNO = "./input/hg38_dfam.nrph.hits"
else:
    # DEV
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"
    DFAM_ANNO = "D:/000_WORK/ParkJiHye/20200914/hg38_dfam.nrph.hits"  # 14
    # DFAM_ANNO = "D:/000_WORK/ParkJiHye/20200914/hg38_dfam.nrph.hits"

PROJECT_NAME = WORK_DIR.split("/")[-2]
FILTERED_CDS_INFO = "filtered_hg38_refFlat.txt"

# name, pam_seq, len_spacer, win_size_arr
INIT = ['TE_trgt', 'NGG', 20, [10, 10]]

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
############### end setting env #################


def split_file_step_0():
    util = Util.Utils()
    util.split_big_file_to_files(DFAM_ANNO, 14, 400000)  # nrph
    # util.split_big_file_to_files(DFAM_ANNO, 296, 400000)  # 309 with unapproved chromosome


def step_1():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    fl_num = 0
    dfam_info = util.read_csv_ignore_N_line(DFAM_ANNO + str(fl_num), '\t', 0)

    # if SYSTEM_NM != 'Linux':
    #     dfam_info = dfam_info[:10]

    dfam_dict = logic_prep.make_list_to_dict_by_ele_as_key(dfam_info, 0)

    # splited_dfam_list = np.array_split(dfam_info, MULTI_CNT)
    #
    # print("platform.system() : ", SYSTEM_NM)
    # print("total cpu_count : ", str(TOTAL_CPU))
    # print("will use : ", str(MULTI_CNT))
    # pool = mp.Pool(processes=MULTI_CNT)
    #
    # pool_list = pool.map(get_trgt, splited_dfam_list)

    header = ['chr', 'tot_seq', 'fam_nm', 'index', 'strand', 'trns_flag']
    for key, val_list in dfam_dict.items():
        result_list = get_trgt(val_list)
        util.make_excel(WORK_DIR + "output/TE_trgt_" + str(fl_num) + "_" + key, header, result_list)


def multi_step_1():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    fl_num = 0
    dfam_info = util.read_csv_ignore_N_line(DFAM_ANNO + str(fl_num), '\t', 0)

    if SYSTEM_NM != 'Linux':
        dfam_info = dfam_info[:100]

    dfam_dict = logic_prep.make_list_to_dict_by_ele_as_key(dfam_info, 0)

    header = ['chr', 'tot_seq', 'fam_nm', 'index', 'strand', 'trns_flag']
    for key, val_list in dfam_dict.items():
        splited_dfam_list = np.array_split(val_list, MULTI_CNT)

        print("platform.system() : ", SYSTEM_NM)
        print("total cpu_count : ", str(TOTAL_CPU))
        print("will use : ", str(MULTI_CNT))
        pool = mp.Pool(processes=MULTI_CNT)

        pool_list = pool.map(get_trgt, splited_dfam_list)
        result_list = logic_prep.merge_multi_list(pool_list)
        print(type(result_list))

        util.make_excel(WORK_DIR + "output/TE_trgt_" + str(fl_num) + "_" + key, header, result_list)

"""
0 #seq_name	
1 family_acc	
2 family_name	
3 bits	
4 e-value	
5 bias	
6 hmm-st	
7 hmm-en	
8 strand	
9 ali-st	
0 ali-en	
1 env-st 11	
2 env-en 12	
3 sq-len	
4 kimura_div
"""
def get_trgt(dfam_list):
    def_nm = "get_trgt"
    print("multi_processing >>>", def_nm)
    util = Util.Utils()
    logic = Logic.Logics()
    result_list = []

    file_nm_arr = ['chrX', 'chrY']
    for f_num in range(1, 23):
        file_nm_arr.append("chr" + str(f_num))

    pam_seq = INIT[1]
    len_spacer = INIT[2]
    len_f_win = INIT[3][0]
    len_b_win = INIT[3][1]
    len_pam = len(pam_seq)

    chr_nm = dfam_list[0][0]
    p_seq, m_seq = util.read_file_by_biopython(REF_DIR + chr_nm + ".fa", "fasta")

    cds_info = util.read_csv_ignore_N_line(WORK_DIR + "input/" + FILTERED_CDS_INFO, "\t")

    for dfam in dfam_list:
        fam_nm = dfam[2]
        strand = dfam[8]

        env_st = 0
        env_en = 0
        if strand == "+":
            env_st += int(dfam[11])
            env_en += int(dfam[12])
        else:
            env_st += int(dfam[12])
            env_en += int(dfam[11])

        trgt_p_seq = p_seq[env_st: env_en]
        trgt_m_seq = m_seq[env_st: env_en]
        # len_trgt = len(trgt_p_seq)

        p_trgt_seq_f = p_seq[env_st - len_f_win: env_st]
        p_trgt_seq_b = p_seq[env_en: env_en + len_b_win]

        m_trgt_seq_f = m_seq[env_en: env_en + len_f_win]
        m_trgt_seq_b = m_seq[env_st - len_b_win: env_st]

        for i in range(len(trgt_p_seq) - len_spacer - len_pam + 1):
            p_pam = trgt_p_seq[i + len_spacer: i + len_spacer + len_pam]
            m_pam = trgt_m_seq[i: i + len_pam]

            real_pos_st = env_st + i
            real_pos_en = env_st + i + len_spacer + len_pam

            # check + strand
            if logic.match(0, p_pam, pam_seq):
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

                trns_flag = False
                for cds_arr in cds_info:
                    # gen_sym = cds_arr[0]
                    # nm_id = cds_arr[1]
                    trns_st = int(cds_arr[4])
                    trns_en = int(cds_arr[5])
                    if trns_st < real_pos_st and real_pos_en < trns_en:
                        trns_flag = True
                        break

                result_list.append([chr_nm, tot_seq, fam_nm, str(real_pos_st) + "-" + str(real_pos_en), '+', trns_flag])

            # check - strand
            if logic.match(0, m_pam, pam_seq[::-1]):
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

                trns_flag = False
                for cds_arr in cds_info:
                    # gen_sym = cds_arr[0]
                    # nm_id = cds_arr[1]
                    trns_st = int(cds_arr[4])
                    trns_en = int(cds_arr[5])
                    if trns_st < real_pos_st and real_pos_en < trns_en:
                        trns_flag = True
                        break

                result_list.append([chr_nm, tot_seq, fam_nm, str(real_pos_st) + "-" + str(real_pos_en), '-', trns_flag])

    print("DONE multi_processing >>>", def_nm)
    return result_list


def test():
    util = Util.Utils()
    cds_info = util.read_csv_ignore_N_line(WORK_DIR + "input/" + FILTERED_CDS_INFO, "\t")
    real_pos_st = 89223770
    real_pos_en = 89223793
    for cds_arr in cds_info:
        if 'chr8' != cds_arr[2]:
            continue

        trns_st = int(cds_arr[4])
        trns_en = int(cds_arr[5])
        if trns_st < real_pos_st and real_pos_en < trns_en:
            print("here")


def test_1000():
    path = "D:/000_WORK/ParkJiHye/20210330/Human Full Genome_TandemRepeat_TRD_20210330.txt"

    with open(path, 'r') as f:
        for i in range(10):
            print(f.readline().replace("\n", ""))


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # split_file_step_0()
    # step_1()
    # multi_step_1()
    test_1000()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
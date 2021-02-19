# TE_targeting_DNA_clock

    [Procedure]
        1. analyze by MultiProcessing_TE_0.py
            input : DFAM_ANNO = "D:/000_WORK/ParkJiHye/20200914/hg38_dfam.nrph.hits"
                    FILTERED_CDS_INFO = "filtered_hg38_refFlat.txt"
                    
            1-1. split DFAM_ANNO file by split_file_step_0()
            1-2. analyze splited DFAM_ANNO files from 1-1 by multi_step_1()
            
        2. merge by MultiProcessing_TE_1_fl1_by_1.py
            2-1.
            
            2-n. make_ACGT_NN()

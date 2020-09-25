
class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"
        self.strt_NM_ = "NM_"

    def sort_list_by_ele(self, data_list, ele_idx, up_down_flag=True):
        result_list = []
        for tmp_arr in sorted(data_list, key=lambda tmp_arr: tmp_arr[ele_idx], reverse=up_down_flag):
            result_list.append(tmp_arr)
        return result_list

    def make_list_to_dict_by_ele_as_key(self, dfam_list, ele_key):
        result_dict = {}
        for dfam_arr in dfam_list:
            tmp_key = dfam_arr[ele_key]
            if tmp_key in result_dict:
                result_dict[tmp_key].append(dfam_arr)
            else:
                result_dict.update({tmp_key: [dfam_arr]})
        return result_dict

    def merge_multi_list(self, pool_list):
        result_list = []
        for split_list in pool_list:
            result_list.extend(split_list)

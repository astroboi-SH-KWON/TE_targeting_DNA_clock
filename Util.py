import glob
import openpyxl

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    """
    get file lists in target dir by target ext
    :param
        path : target dir + "*." + target ext
    :return
        ['target dir/file_name.target ext', 'target dir/file_name.target ext' ...]
    """
    def get_files_from_dir(self, path):
        return glob.glob(path)

    def read_csv_ignore_N_line(self, path, deli_str=",", n_line=1):
        result_list = []
        with open(path, "r") as f:
            for ignr_line in range(n_line):
                header = f.readline()
                print(header)
            while True:
                tmp_line = f.readline().replace("\n", "")
                if tmp_line == '':
                    break

                result_list.append(tmp_line.split(deli_str))
        return result_list

    def make_excel_row(self, sheet, row, data_arr, col=1):
        for idx in range(len(data_arr)):
            sheet.cell(row=row, column=(col + idx), value=data_arr[idx])

    def make_excel(self, path, header, data_list, strt_idx=0):
        workbook = openpyxl.Workbook()
        sheet = workbook.active

        row = 1
        self.make_excel_row(sheet, row, header[strt_idx:])

        for data_arr in data_list:
            row += 1
            self.make_excel_row(sheet, row, data_arr[strt_idx:])

        workbook.save(filename=path + self.ext_xlsx)

    def make_csv(self, path, header, data_list, strt_idx=0, deli=','):
        with open(path, 'a') as f:
            tmp_head = ''
            for head in header[strt_idx:]:
                tmp_head += (head + deli)
            f.write(tmp_head[:-1] + "\n")

            for data_arr in data_list:
                tmp_row = ''
                for row_val in data_arr[strt_idx:]:
                    tmp_row += (str(row_val) + deli)
                f.write(tmp_row[:-1] + "\n")
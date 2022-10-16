# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 10:58:41 2022

@author: tianb
"""

'''this code is used to delete data files under the 'duc_result_data_file' directory on hpc platform
and all files in folders under it. We define a recursive function to do this'''

import shutil
import os
from pathlib import Path

def del_files(dir_path):
    if os.path.isfile(dir_path):
        try:
            os.remove(dir_path)
        except BaseException as e:
            print(e)
    elif os.path.isdir(dir_path):
        file_list = os.listdir(dir_path)
        for file_name in file_list:
            this_file = os.path.join(dir_path, file_name)
            del_files(this_file)

#dir_path = 'C:\\Apps\\duc_result_data_file'
dir_path = '/hpc/puhome/21040136r/my_programs/DUC_Programs/duc_result_data_file'

if __name__ == '__main__':
    del_files(dir_path)
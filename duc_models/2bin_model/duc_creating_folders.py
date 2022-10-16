# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 10:01:29 2022

@author: tianb
"""

'''This code is used for making new folders for storing data result on hpc platform'''

import os

#set the path where you want to create a new folder
path = '/hpc/puhome/21040136r/my_programs/DUC_Programs/test_result_data_file/'

#create 20 new folders for storing the data result of 20 instances
for i in range(1, 21):
    os.makedirs(path + 'instance{0}'.format(i), exist_ok = True)
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:07:59 2022

@author: tianb
"""

'''
This code is for solving deterministic Unit-Commitment problem with ramp 
up/down constraints. Strong valid inequalities are used to speed up the algorithm.
'''

import numpy as np
import pandas as pd
import sys
from duc_model_funcs import *
from duc_params import *
       

if __name__ == '__main__':
    '''there are totally 20 instances and we run each instance with four algorithm and three times for each'''
    start = int(sys.argv[1])
    flag = int(sys.argv[2])
    for current_instance_index in range(start, start+1):
        #for each instance we run four kinds of algorithm: original, with only user cuts to cutpool, with only
        #callback, with both cutpool and callback. And for each kind we run 3 times
        for current_num_of_run in range(1, 4):
            if flag == 0:
                #print('without any cuts:')
                solve_model(usercutpool = False, usercutcallback = False, instance_index = current_instance_index, \
                        num_of_run = current_num_of_run)
            elif flag == 1:
                #print('with only usercutpool cuts:')
                solve_model(usercutpool = True, usercutcallback = False, instance_index = current_instance_index, \
                        num_of_run = current_num_of_run)
            elif flag == 2:
                #print('with only callback:')
                solve_model(usercutpool = False, usercutcallback = True, instance_index = current_instance_index, \
                        num_of_run = current_num_of_run)
            else:
                #print('with both usercutpool and callback:')
                solve_model(usercutpool = True, usercutcallback = True, instance_index = current_instance_index, \
                        num_of_run = current_num_of_run)
    #set data result file path
    uc_data_result_file = filepath + '/uc_data_result_{0}_{1}'.format(start, flag) + '.csv'
    #export the file
    export_timegap_info(uc_data_result_file)
        
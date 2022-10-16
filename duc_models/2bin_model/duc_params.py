# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 21:54:03 2022

@author: tianb
"""
import numpy as np
import pandas as pd

''''initialization of all parameters used in the model'''
#reading default parameters
file = 'duc_parameters.csv'
data = pd.read_csv(file, header = 0, index_col = 0)
#operation period, here we set it to be 24 hours
uc_opt_period = np.array(data.iloc[0, 0], dtype = int)

#number of generators for 8 types of generator, instance1 starts from row 17 and instance 20 end in row 37
#20*8 array
uc_gen_num = np.array(data.iloc[18:38, 0:8], dtype = int)
#calculating the total number of generators for all 20 instances, 1*20 array
uc_total_gen_num = uc_gen_num.sum(axis = 1)
#minimum on time period when a generator starts up, for 8 types of generators
uc_min_on = np.array(data.iloc[1, 0:8], dtype = int)
#minimum off time period when a generator starts up, for 8 types of generators
uc_min_off = np.array(data.iloc[2, 0:8], dtype = int)
#generation upper bound for 8 types of generators
uc_gen_ub = np.array(data.iloc[3, 0:8], dtype = float)
#generation lower bound for 8 types of generators
uc_gen_lb = np.array(data.iloc[4, 0:8], dtype = float)
#star up upper bound for 8 types of generators
uc_start_ub = np.array(data.iloc[5, 0:8], dtype = float)
#ramp up upper bound for 8 types of generators
uc_rampup_ub = np.array(data.iloc[6, 0:8], dtype = float)
#ramp down upper bound for 8 types of generators
uc_rampdown_ub = np.array(data.iloc[7, 0:8], dtype = float)
#unit start up cost ($/h) for 8 types of generators 
uc_unit_startup_cost = np.array(data.iloc[8, 0:8], dtype = float)
#unit shut down cost ($/h) for 8 types of generators
uc_unit_shutdown_cost = np.array(data.iloc[9, 0:8], dtype = float)
#number of break points we use to approximate the quadratic generation cost f(x_t) by a piecewise linear function
uc_gen_cost_pwl_num = np.array(data.iloc[10, 0:8], dtype = int) 
#the coefficient of the second order term in the quadratic generation cost function f(x_t)
uc_gen_cost_func_a = np.array(data.iloc[11, 0:8], dtype = float)
#the coefficient of the first order term in the quadratic generation cost function f(x_t)
uc_gen_cost_func_b = np.array(data.iloc[12, 0:8], dtype = float)
#the value of the constant term in the quadratic generation cost function f(x_t)
uc_gen_cost_func_c = np.array(data.iloc[13, 0:8], dtype = float)  
#spinning reserve percentage
uc_spin_reserve = np.array(data.iloc[15, 0:24], dtype = float)

#set 24 periods load for each insatnce , 20*24 array
uc_load = np.zeros((20,24), dtype = float)
for instance_index in range(0,20):
    #total system load capacity (total demand capacity) for 20 instances
    uc_total_load_capacity = sum(uc_gen_ub[i]*uc_gen_num[instance_index, i] for i in range(0,8))
    #load (electricity demand) for 24 time period
    uc_load[instance_index, 0:24] = uc_total_load_capacity*np.array(data.iloc[14, 0:24]) 

#calculate the accumulation number of first k types of generators for 20 instances, because we need to use the 
#start and end position of each kind of generators in [1, G]
#uc_gen_accumu_num is the accumulation number of 8 kinds of generators for each instance, it is a 20*9 array because
#the first column is 0
uc_gen_accumu_num = np.zeros((20,9), dtype = int)
for instance_index in range(0,20):
    start_value = 0
    for k in range(1,9):
        start_value += int(uc_gen_num[instance_index, k-1])
        uc_gen_accumu_num[instance_index, k] = start_value

#time limit on solving
time_limit = np.array(data.iloc[16, 0])
#MIP gap tolerance on solving
mip_gap_tolerance = np.array(data.iloc[17, 1])
#set file output path
#filepath = 'C:\\Apps\\duc_result_data_file'
filepath = '/hpc/puhome/21040136r/my_programs/DUC_Programs/test_result_data_file'
#filepath = '/hpc/puhome/21040136r/my_programs/DUC_Programs/test_data_file'

#set data result array, after all optimization we write this array to uc_data_result.csv in the main func
uc_data_result = np.zeros((1,9), dtype = float)
#count the number of usercut callback added
uc_num_usercut_callback_added = 0
#calculating error tolerance
uc_eps = 10e-5 

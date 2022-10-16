# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 10:04:35 2022

@author: tianb
"""

'''
author: Bin Tian
date: 2/8/2022
'''
import numpy as np
import pandas as pd
from duc_params import *

#array used to store any violation index of g, x_t, and x_[t-k]
a7_violation = np.array([0,0,0])
a7_lhs_rhs = np.array([0.0,0.0], dtype = float)
b7_violation = np.array([0,0,0])
b7_lhs_rhs = np.array([0.0,0.0], dtype = float)
a8_violation = np.array([0,0,0])
a8_lhs_rhs = np.array([0.0,0.0], dtype = float)
b8_violation = np.array([0,0,0])
b8_lhs_rhs = np.array([0.0,0.0], dtype= float)
eps = 10e-5
'''checking violation'''
def violation_check(x_var, y_var, instance_index):
    start_index = 0
    end_index = 0
    #if usercutpool we add usercuts to the usercutpool
    #add usercuts for 8 types of generators one by one
    for gen_type in range(0, 8):
        #if the number of generators of current type is not 0
        if uc_gen_num[instance_index-1, gen_type] > 0:
            #update the end_index to the end of current type of generators
            end_index += uc_gen_num[instance_index-1, gen_type]
            #constraints(7a) & (7b): ramp up and down constraint for two vairbales
            for k in range(1, uc_opt_period-1):
                if uc_gen_ub[gen_type] - uc_gen_lb[gen_type] - k*uc_rampup_ub[gen_type] > 0:
                    #constriant(7a): ramp up constraint for two variables
                    #in the names, dvar for double var
                    for g in range(start_index+1, end_index+1):
                        for t in range(k+1, uc_opt_period):
                            lhs = x_var['x_{0}_{1}'.format(g,t)] - x_var['x_{0}_{1}'.format(g,t-k)]
                            rhs = (uc_gen_lb[gen_type] + k*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])* \
                                y_var['y_{0}_{1}'.format(g,t+1)] + uc_start_ub[gen_type]*y_var['y_{0}_{1}'.format(g, t)] \
                                - uc_gen_lb[gen_type]*y_var['y_{0}_{1}'.format(g, t-k)]
                            if lhs > rhs:
                                global a7_lhs_rhs
                                a7_lhs_rhs = np.vstack([a7_lhs_rhs, [lhs, rhs]])
                                global a7_violation
                                a7_violation = np.vstack([a7_violation, [g, t, t-k]])
                        for t in range(2, uc_opt_period-k+1):
                            lhs = x_var['x_{0}_{1}'.format(g,t)] - x_var['x_{0}_{1}'.format(g,t+k)]
                            rhs = (uc_gen_lb[gen_type] + k*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])* \
                                y_var['y_{0}_{1}'.format(g,t-1)] + uc_start_ub[gen_type]*y_var['y_{0}_{1}'.format(g, t)] \
                                - uc_gen_lb[gen_type]*y_var['y_{0}_{1}'.format(g, t+k)]
                            if lhs > rhs+eps:
                                global b7_lhs_rhs
                                b7_lhs_rhs = np.vstack([b7_lhs_rhs, [lhs, rhs]])
                                global b7_violation
                                b7_violation = np.vstack([b7_violation, [g, t, t-k]])
            for k in range(1, uc_opt_period):
                if uc_gen_ub[gen_type] - uc_gen_lb[gen_type] - k*uc_rampup_ub[gen_type] > 0:
                    #constraint(8a): ramp up constraint for two variables
                    for g in range(start_index+1, end_index+1):
                        for t in range(k+1, uc_opt_period+1):
                            lhs = x_var['x_{0}_{1}'.format(g,t)] - x_var['x_{0}_{1}'.format(g,t-k)]
                            rhs = (uc_gen_lb[gen_type] + k*uc_rampup_ub[gen_type])*y_var['y_{0}_{1}'.format(g,t)] -\
                                uc_gen_lb[gen_type]*y_var['y_{0}_{1}'.format(g, t-k)]
                            if lhs > rhs+eps:
                                global a8_lhs_rhs
                                a8_lhs_rhs = np.vstack([a8_lhs_rhs, [lhs, rhs]])
                                global a8_violation
                                a8_violation = np.vstack([a8_violation, [g, t, t-k]])
                        for t in range(1, uc_opt_period-k+1):
                            lhs = x_var['x_{0}_{1}'.format(g,t)] - x_var['x_{0}_{1}'.format(g,t+k)]
                            rhs = (uc_gen_lb[gen_type] + k*uc_rampup_ub[gen_type])*y_var['y_{0}_{1}'.format(g,t)] -\
                                uc_gen_lb[gen_type]*y_var['y_{0}_{1}'.format(g, t+k)]
                            if lhs > rhs+eps:
                                global b8_lhs_rhs
                                b8_lhs_rhs = np.vstack([b8_lhs_rhs, [lhs, rhs]])
                                global b8_violation
                                b8_violation = np.vstack([b8_violation, [g, t, t-k]])
            
            start_index = end_index

'''instance'''
instance_index = 1
#calculating the number of decision variables x_var
var_num = uc_total_gen_num[instance_index-1]*uc_opt_period
'''reading solution'''
x_var = {}
y_var = {}
#file = 'MILP_solution_instance{0}_1.csv'.format(instance_index)
file = 'Usercuts_solution_instance{0}_1.csv'.format(instance_index)
#file = 'Usercuts_Callback_solution_instance{0}_1.csv'.format(instance_index)
#file = 'Usercuts_Callback_solution_instance{0}_1.csv'.format(instance_idnex)
data = pd.read_csv(file, header = 0, index_col = None, nrows = var_num*2)
for i in range(0, var_num):
    key, value = data.iloc[i,:]
    y_var[key] = value
for i in range(var_num, 2*var_num):
    key, value = data.iloc[i,:]
    x_var[key] = value

if __name__ == '__main__':
    
    violation_check(x_var, y_var, instance_index)
    
    if len(a7_lhs_rhs) > 2:
        print('the violation of constriant (7a):')
        for i in range(1, len(a7_lhs_rhs)):
            print('lhs: {0}'.format(a7_lhs_rhs[i,0]), '<=', 'rhs: {0}'.format(a7_lhs_rhs[i,1]))
    else:
        print('no violation of constraint (7a)')
        
    if len(b7_lhs_rhs) > 2:
        print('the violation of constriant (7b):')
        for i in range(1, len(b7_lhs_rhs)):
            print('lhs: {0}'.format(b7_lhs_rhs[i,0]), '<=', 'rhs: {0}'.format(b7_lhs_rhs[i,1]))
    else:
        print('no violation of constraint (7b)')
        
    if len(a8_lhs_rhs) > 2:
        print('the violation of constriant (8a):')
        for i in range(1, len(a8_lhs_rhs)):
            print('lhs: {0}'.format(a8_lhs_rhs[i,0]), '<=', 'rhs: {0}'.format(a8_lhs_rhs[i,1]))
    else:
        print('no violation of constraint (8a)')
        
    if len(b8_lhs_rhs) > 2:
        print('the violation of constriant (8b):')
        for i in range(1, len(b8_lhs_rhs)):
            print('lhs: {0}'.format(b8_lhs_rhs[i,0]), '<=', 'rhs: {0}'.format(b8_lhs_rhs[i,1]))
    else:
        print('no violation of constraint (8b)')
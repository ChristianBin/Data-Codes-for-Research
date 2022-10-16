# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 21:23:49 2022

@author: tianb
"""
import numpy as np
import pandas as pd
from duc_params import *

def set_file_name(usercutpool, usercutcallback):
    '''func to set file name according to what optimization strategy you are using'''
    
    #set the logfile path and name, write the logfile to csv file
    if not usercutpool and not usercutcallback:
        filename = 'MILP'
    elif usercutpool and not usercutcallback:
        filename = 'Usercuts'
    elif not usercutpool and usercutcallback:
        filename = 'Callback'
    else:
        filename = 'Usercuts_Callback'
    return filename

def export_sol_result(duc_model, usercutpool, usercutcallback, instance_index, num_of_run): 
    '''export the solution result and store the data_result'''
    
    #get the data results that we want
    start_shut_cost = duc_model.kpi_value_by_name('total start up and shut down cost')
    generation_cost = duc_model.kpi_value_by_name('total genration cost')
    obj_value = duc_model.objective_value
    cpu_time = duc_model.solve_details.time
    mip_gap = duc_model.solve_details.gap
    total_usercuts_num = duc_model.number_of_user_cut_constraints
    num_nodes_explored = duc_model.solve_details.nb_nodes_processed
    usercuts_applied = duc_model.cplex.solution.MIP.get_num_cuts(duc_model.cplex.solution.MIP.cut_type.user)
    best_bound = duc_model.solve_details.best_bound

    #set the solution and model file path and name, write the solution ot csv file, export the model as 
    #.lp file
    sol_filename = set_file_name(usercutpool, usercutcallback) + '_solution'
    sol_filename = filepath + '/instance{0}/'.format(instance_index) + sol_filename \
        + '_instance{0}_{1}.csv'.format(instance_index, num_of_run)
    #we store the solution to an array and the write to a csv file
    sol_name = np.array([])
    sol_value = np.array([])
    for v in duc_model.iter_binary_vars():
        sol_name = np.append(sol_name, str  (v))
        sol_value = np.append(sol_value, str(v.solution_value))
    for v in duc_model.iter_continuous_vars():
        sol_name = np.append(sol_name, str(v))
        sol_value = np.append(sol_value, str(v.solution_value))
    sol_name_value = np.vstack([sol_name, sol_value])
    sol_name_value = pd.DataFrame(sol_name_value.T)
    sol_name_value.to_csv(sol_filename, index = False, header = ['variable', 'value'])
            
    #store the data result to uc_data_result
    global uc_data_result
    uc_data_result = np.vstack([uc_data_result, [start_shut_cost, generation_cost, obj_value, best_bound, cpu_time, \
        mip_gap, total_usercuts_num, usercuts_applied, num_nodes_explored]])
        
def export_timegap_info(uc_data_result_file):
    '''export the final data_result'''
    
    global uc_data_result
    uc_data_result = pd.DataFrame(uc_data_result)
    uc_data_result_header = ['total_start_shut_cost','total_gen_cost', 'obj_value', 'best_bound', 'cpu_time(s)', \
                             'mip_gap', 'total_usercuts_num', 'total_usercuts_applied_num', \
                                'nodes_explored_num']  
    uc_data_result.to_csv(uc_data_result_file, index = True, header = uc_data_result_header)
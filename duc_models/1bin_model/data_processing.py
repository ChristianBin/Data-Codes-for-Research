# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 10:44:58 2022

@author: tianb
"""

'''This code is for processing the resulting data of DUC model'''
import pandas as pd
import numpy as np


def file_processing(file_name):
    data = pd.read_csv(file_name, header = 0, index_col = 0, skiprows = [0])
    data_column_num = data.shape[1]
    column_mean = np.array([data.iloc[:,i].mean() for i in range(0, data_column_num)])
    return column_mean
    
#set the path
file_path = 'C:\\Bin Files\The Hong Kong Polytechnic University\\Research\\Unit-Commitment\\working papers\codes\\duc_python_codes\\Exp1-8_type_generators\\data_result\\'
#set the path for specific data set
file_path += '2bin_model_result\\result2\\gen2_32_threads\\1bin_user_results\\add_with_sq\\'
#use to store final data
data_result = np.array([])
#set the index for final result data file
data_result_index = np.array([])
#set the header for final result data file
data_result_header = ['instances', 'total_start_shut_cost','total_gen_cost', 'obj_value', 'best_bound', 'cpu_time(s)', \
                      'mip_gap', 'total_usercuts_num', 'total_usercuts_applied_num', 'nodes_explored_num']
#set the name for final result
final_data_result = file_path +'20_instance_data.csv'
    
if __name__ == '__main__':
    #for 20 instances
    for i in range(1, 11):
        #for 4 algorithms
        for j in range(1, 2):
            file_name = file_path + 'uc_data_result_{0}_{1}.csv'.format(i,j)
            new_row = file_processing(file_name)
            if len(data_result) == 0:
                data_result = np.append(data_result, new_row)
            else:
                data_result = np.vstack([data_result, new_row])
            #append the row index
            data_result_index = np.append(data_result_index,'instance_{0}_{1}'.format(i,j))
    
    #export the file
    data_result_index = np.array([data_result_index])
    data_result = np.concatenate((data_result_index.T, data_result), axis=1)
    data_result = pd.DataFrame(data_result)
    data_result.to_csv(final_data_result, index = False, header = data_result_header)
    
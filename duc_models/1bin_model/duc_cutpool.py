# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 10:40:01 2022

@author: tianb
"""
from docplex.mp.model import Model
from duc_params import *
import math
    
def add_cutpool_cuts(model, x_var, y_var, instance_index):
    #use start_index and end_index to denote the accumulating number of each type of generators
    #for current instance
    start_index = 0
    end_index = 0
    #if usercutpool we add usercuts to the usercutpool
    #add usercuts for 8 types of generators one by one
    for gen_type in range(0, 8):
        #if the number of generators of current type is not 0
        if uc_gen_num[instance_index-1, gen_type] > 0:
            #update the end_index to the end of current type of generators
            end_index += uc_gen_num[instance_index-1, gen_type]
            
            '''add 1-binary usercuts to the 2bin model'''
            
            #the biggest number of ramping from V_bar to C_bar for current type of generators
            ramp_nb_limit = math.floor((uc_gen_ub[gen_type] - uc_start_ub[gen_type])/uc_rampup_ub[gen_type])
            
            #constraint (3a) & (3b)
            k_prime = min(ramp_nb_limit, uc_opt_period-2)
            for k in range(0, k_prime+1):
                #constraint (3a): rampup upper bound constraint
                model.add_user_cut_constraints((x_var[g, t] <= uc_gen_ub[gen_type]*y_var[g, t] - sum((uc_gen_ub[gen_type]-\
                    uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])*(y_var[g, t-i] - y_var[g, t-i-1]) \
                    for i in range(0, k+1)) for g in range(start_index+1, end_index+1) \
                    for t in range(k+2, uc_opt_period+1)), names = '1bin_usercut_3a')
            
                #constraint (3b): rampdown upper bound constraint
                model.add_user_cut_constraints((x_var[g, t] <= uc_gen_ub[gen_type]*y_var[g, t] - sum((uc_gen_ub[gen_type]-\
                    uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])*(y_var[g, t+i] - y_var[g, t+i+1]) \
                    for i in range(0, k+1)) for g in range(start_index+1, end_index+1) \
                    for t in range(1, uc_opt_period-k)), names = '1bin_usercut_3b')
                
                        
            #usercut constraint (5a) & (5b)
            beta = min(uc_min_on[gen_type]-1, (uc_gen_ub[gen_type] - uc_start_ub[gen_type])/uc_rampup_ub[gen_type])
            gamma = min(ramp_nb_limit, uc_opt_period-3)
            for k in range(0, gamma+1):
                #constraint(5a): rampup upper bound constraint
                model.add_user_cut_constraints((x_var[g,t] <= (uc_gen_ub[gen_type] - beta*uc_rampup_ub[gen_type])*y_var[g,t] +\
                    beta*uc_rampup_ub[gen_type]*y_var[g, t+1] - sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] -\
                    i*uc_rampup_ub[gen_type])*(y_var[g, t-i] - y_var[g, t-i-1]) for i in range(0, k+1)) \
                    for g in range(start_index+1, end_index+1) for t in range(k+2, uc_opt_period)), names = '1bin_usercut_5a')
                 
                #constraint(5b): rampdown upper bound constraint
                model.add_user_cut_constraints((x_var[g,t] <= (uc_gen_ub[gen_type] - beta*uc_rampup_ub[gen_type])*y_var[g,t] +\
                    beta*uc_rampup_ub[gen_type]*y_var[g, t-1] - sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] -\
                    i*uc_rampup_ub[gen_type])*(y_var[g, t+i] - y_var[g, t+i+1]) for i in range(0, k+1)) \
                    for g in range(start_index+1, end_index+1) for t in range(2, uc_opt_period-k)), names = '1bin_usercut_5b')
                   
            #usercut constraint (6a) & (6b)
            beta1 = 0
            beta2 = min(uc_min_on[gen_type], (uc_gen_ub[gen_type] - uc_start_ub[gen_type])/uc_rampup_ub[gen_type])
            gamma = min(ramp_nb_limit, uc_opt_period-2)
            for k in range(1, gamma+1):
                #constraint(6a): rampup upper bound constraint
                #add usercut (6a) with beta = 0
                model.add_user_cut_constraints((x_var[g, t] <= (uc_start_ub[gen_type] + beta1*uc_rampup_ub[gen_type])*y_var[g, t] +\
                    (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - beta1*uc_rampup_ub[gen_type])*y_var[g, t-1] -\
                    sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])*(y_var[g, t-i] -\
                    y_var[g, t-i-1]) for i in range(1, k+1)) for g in range(start_index+1, end_index+1) \
                    for t in range(k+2, uc_opt_period+1)), names = '1bin_usercut_6a_I')
                #add usercut (6a) with beta = beta2
                model.add_user_cut_constraints((x_var[g, t] <= (uc_start_ub[gen_type] + beta2*uc_rampup_ub[gen_type])*y_var[g, t] +\
                    (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - beta2*uc_rampup_ub[gen_type])*y_var[g, t-1] -\
                    sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])*(y_var[g, t-i] -\
                    y_var[g, t-i-1]) for i in range(1, k+1)) for g in range(start_index+1, end_index+1) \
                    for t in range(k+2, uc_opt_period+1)), names = '1bin_usercut_6a_II')
                    
                ##constraint(6b): rampdown upper bound constraint
                #add usercut (6b) with beta = 0
                model.add_user_cut_constraints((x_var[g, t] <= (uc_start_ub[gen_type] + beta1*uc_rampup_ub[gen_type])*y_var[g, t] +\
                    (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - beta1*uc_rampup_ub[gen_type])*y_var[g, t+1] -\
                    sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])*(y_var[g, t+i] -\
                    y_var[g, t+i+1]) for i in range(1, k+1)) for g in range(start_index+1, end_index+1) \
                    for t in range(1, uc_opt_period-k)), names = '1bin_usercut_6b_I')
                #add usercut (6b) with beta = beta2
                model.add_user_cut_constraints((x_var[g, t] <= (uc_start_ub[gen_type] + beta2*uc_rampup_ub[gen_type])*y_var[g, t] +\
                    (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - beta2*uc_rampup_ub[gen_type])*y_var[g, t+1] -\
                    sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])*(y_var[g, t+i] -\
                    y_var[g, t+i+1]) for i in range(1, k+1)) for g in range(start_index+1, end_index+1) \
                    for t in range(1, uc_opt_period-k)), names = '1bin_usercut_6b_II')
            
            #constraints(7a) & (7b): ramp up and down constraint for two vairbales
            for k in range(1, uc_opt_period-1):
                if uc_gen_ub[gen_type] - uc_gen_lb[gen_type] - k*uc_rampup_ub[gen_type] > 0:
                    #constriant(7a): ramp up constraint for two variables
                    #in the names, dvar for double var
                    
                    model.add_user_cut_constraints((x_var[g,t] - x_var[g,t-k] <=\
                        (uc_gen_lb[gen_type] + k*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*y_var[g,t+1] + \
                        uc_start_ub[gen_type]*y_var[g,t] - uc_gen_lb[gen_type]*y_var[g,t-k] for g in \
                        range(start_index+1, end_index+1) for t in range(k+1, uc_opt_period)), \
                        names = 'usercut_rampup_seven_a')
                      
                    #constraint(7b): ramp down constraint for two variables
                    model.add_user_cut_constraints((x_var[g,t] - x_var[g,t+k] <=\
                        (uc_gen_lb[gen_type] + k*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])*y_var[g,t-1] + \
                        uc_start_ub[gen_type]*y_var[g,t] - uc_gen_lb[gen_type]*y_var[g,t+k] for g in \
                        range(start_index+1, end_index+1) for t in range(2, uc_opt_period-k+1)), \
                        names = 'usercut_rampdown_seven_b')
            
            #constraints (8a) & (8b): ramp up and down constraint for two variables
            for k in range(1, uc_opt_period):
                if uc_gen_ub[gen_type] - uc_gen_lb[gen_type] - k*uc_rampup_ub[gen_type] > 0:
                    
                    #constraint(8a): ramp up constraint for two variables
                    model.add_user_cut_constraints((x_var[g,t] - x_var[g,t-k] <= (uc_gen_lb[gen_type] + \
                        k*uc_rampup_ub[gen_type])*y_var[g,t] - uc_gen_lb[gen_type]*y_var[g,t-k] for g in \
                        range(start_index+1, end_index+1) for t in range(k+1, uc_opt_period+1)), \
                        names = 'usercut_rampup_eight_a')
                        
                    #constraint(8b): ramp down constraint for two variables
                    model.add_user_cut_constraints((x_var[g,t] - x_var[g,t+k] <= (uc_gen_lb[gen_type] + \
                        k*uc_rampdown_ub[gen_type])*y_var[g,t] - uc_gen_lb[gen_type]*y_var[g,t+k] for g in \
                        range(start_index+1, end_index+1) for t in range(1, uc_opt_period-k+1)), \
                        names = 'usercut_rampdown_eight_b')
            
            #constraint (9a) & (9b)
            k_prime = min(uc_opt_period-1, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampup_ub[gen_type]))
            for k in range(1, k_prime+1):
                gamma = min(k-1, uc_min_on[gen_type]-1)
                for sq in range(0, gamma+1):
                    #constraint (9a): rampup upper bound constraint
                    model.add_user_cut_constraints((x_var[g,t] - x_var[g, t-k] <= (uc_gen_lb[gen_type]+k*uc_rampup_ub[gen_type])*\
                        y_var[g,t] - uc_gen_lb[gen_type]*y_var[g, t-k] - sum((uc_gen_lb[gen_type]+(k-i)*uc_rampup_ub[gen_type] -\
                        uc_start_ub[gen_type])*(y_var[g,t-i] - y_var[g, t-i-1]) for i in range(0, sq+1)) \
                        for g in range(start_index+1, end_index+1) for t in range(k+1, uc_opt_period+1)), names = '1bin_usercut_9a')
                    
                    #constraint(9b): rampdown upper bound constraint
                    model.add_user_cut_constraints((x_var[g,t] - x_var[g, t+k] <= (uc_gen_lb[gen_type]+k*uc_rampup_ub[gen_type])*\
                        y_var[g,t] - uc_gen_lb[gen_type]*y_var[g, t+k] - sum((uc_gen_lb[gen_type]+(k-i)*uc_rampup_ub[gen_type] -\
                        uc_start_ub[gen_type])*(y_var[g,t+i] - y_var[g, t+i+1]) for i in range(0, sq+1)) \
                        for g in range(start_index+1, end_index+1) for t in range(1, uc_opt_period-k+1)), names = '1bin_usercut_9b')
                    
            #constraints(10a) & (10b): ramp up and down constraints for two variables
            ub = math.floor((uc_opt_period-1)/2)
            for k in range(1, ub+1):
                if uc_gen_ub[gen_type] - uc_gen_lb[gen_type] - k*uc_rampup_ub[gen_type] > 0:
                    #constraint(10a): ramp up constraint for two variables
                    model.add_user_cut_constraints((x_var[g,t] - x_var[g, t-k] <= uc_start_ub[gen_type]*y_var[g,t] \
                        + uc_rampup_ub[gen_type]*model.sum_vars(y_var[g,t+i] for i in range(1, k)) + \
                        (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*y_var[g, t+k] - \
                        uc_gen_lb[gen_type]*y_var[g,t-k] for g in range(start_index+1, end_index+1) for t in \
                        range(k+1, uc_opt_period-k+1)), names = 'usercut_rmapup_ten_a')
                        
                    #constraint(10b): ramp down constraint for two variables
                    model.add_user_cut_constraints((x_var[g,t] - x_var[g, t+k] <= uc_start_ub[gen_type]*y_var[g,t] \
                        + uc_rampup_ub[gen_type]*model.sum_vars(y_var[g,t-i] for i in range(1, k)) + \
                        (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*y_var[g, t-k] - \
                        uc_gen_lb[gen_type]*y_var[g,t+k] for g in range(start_index+1, end_index+1) for t in \
                        range(k+1, uc_opt_period-k+1)), names = 'usercut_rmapup_ten_b')
                        
            #constraint (11a) & (11b)
            if uc_min_on[gen_type] > 1:
                k_prime = min(uc_opt_period-1, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampup_ub[gen_type]))
                for k in range(1, k_prime+1):
                    gamma = min(k-1, uc_min_on[gen_type]-2)
                    for m in range(0, gamma+1):
                        #for sq in range(0, gamma-m+1):
                        #constraint (11a): rampup upper bound for two variables
                        model.add_user_cut_constraints((x_var[g, t] - x_var[g, t-k] <= (uc_gen_lb[gen_type] +\
                            (k-m)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*y_var[g, t+m+1] +\
                            sum(uc_rampup_ub[gen_type]*y_var[g, t+j] for j in range(1, m+1)) + uc_start_ub[gen_type]*y_var[g,t]-\
                            uc_gen_lb[gen_type]*y_var[g, t-k] - sum((uc_gen_lb[gen_type] + (k-i)*uc_rampup_ub[gen_type] -\
                            uc_start_ub[gen_type])*(y_var[g,t-i] - y_var[g, t-i-1]) for i in range(0, gamma-m+1)) \
                            for g in range(start_index+1, end_index+1) for t in range(k+1, uc_opt_period-m)), names = '1bin_usercut_11a')
                            
                        #constraint (11b): rampdown upper bound for two variables
                        model.add_user_cut_constraints((x_var[g, t] - x_var[g, t+k] <= (uc_gen_lb[gen_type] +\
                            (k-m)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*y_var[g, t-m-1] +\
                            sum(uc_rampup_ub[gen_type]*y_var[g, t-j] for j in range(1, m+1)) + uc_start_ub[gen_type]*y_var[g,t]-\
                            uc_gen_lb[gen_type]*y_var[g, t+k] - sum((uc_gen_lb[gen_type] + (k-i)*uc_rampup_ub[gen_type] -\
                            uc_start_ub[gen_type])*(y_var[g,t+i] - y_var[g, t+i+1]) for i in range(0, gamma-m+1)) \
                            for g in range(start_index+1, end_index+1) for t in range(m+2, uc_opt_period-k+1)), names = '1bin_usercut_11b')
                            
            #constraint (12a) & (12b)
            k_prime = min(uc_opt_period-1, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampup_ub[gen_type]))
            for k in range(1, k_prime+1):
                gamma = min(k-1, uc_min_on[gen_type]-1)
                for m in range(0, gamma+1):
                    #for sq in range(0, gamma-m+1):
                    #constraint(12a): rampup upper bound constraint
                    model.add_user_cut_constraints((x_var[g,t] - x_var[g,t-k] <= sum(uc_rampup_ub[gen_type]*y_var[g,t+j] \
                        for j in range(1,m+1)) + (uc_gen_lb[gen_type]+(k-m)*uc_rampup_ub[gen_type])*y_var[g,t] -\
                        uc_gen_lb[gen_type]*y_var[g,t-k] - sum((uc_gen_lb[gen_type] + (k-i)*uc_rampup_ub[gen_type] -\
                        uc_start_ub[gen_type])*(y_var[g,t-i] - y_var[g,t-i-1]) for i in range(0, gamma-m+1)) \
                        for g in range(start_index+1, end_index+1) for t in range(k+1, uc_opt_period-m+1)), names = '1bin_usercut_12a')
                    
                    #constraint (12b): rampdown upper bound constraint
                    model.add_user_cut_constraints((x_var[g,t] - x_var[g,t+k] <= sum(uc_rampup_ub[gen_type]*y_var[g,t-j] \
                        for j in range(1,m+1)) + (uc_gen_lb[gen_type]+(k-m)*uc_rampup_ub[gen_type])*y_var[g,t] -\
                        uc_gen_lb[gen_type]*y_var[g,t+k] - sum((uc_gen_lb[gen_type] + (k-i)*uc_rampup_ub[gen_type] -\
                        uc_start_ub[gen_type])*(y_var[g,t+i] - y_var[g,t+i+1]) for i in range(0, gamma-m+1)) \
                        for g in range(start_index+1, end_index+1) for t in range(m+1, uc_opt_period-k+1)), names = '1bin_usercut_12b')
            
            start_index = end_index
            
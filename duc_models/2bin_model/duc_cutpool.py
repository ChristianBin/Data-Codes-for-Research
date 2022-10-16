# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 10:40:01 2022

@author: tianb
"""
from docplex.mp.model import Model
from duc_params import *
import math
    
def add_cutpool_cuts(model, x_var, y_var, u_var, instance_index):
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
            
            '''add the 2bin 3period convex hull cust with L=l=2 for L >= 2 of 2bin model to 2bin model as usercuts'''
            '''
            if uc_min_on[gen_type] >= 2:
                #constraint (4): rampup upper bound constraint
                model.add_user_cut_constraints((x_var[g,t] <= uc_start_ub[gen_type]*y_var[g,t] + uc_rampup_ub[gen_type]*(y_var[g,t+1] - u_var[g,t+1]) +\
                    (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - uc_rampup_ub[gen_type])*(y_var[g,t+2] - u_var[g,t+2] - u_var[g,t+1]) \
                    for g in range(start_index+1, end_index+1) for t in range(1, uc_opt_period-1)), names = 'three_period_conv_4')
                
                #constraint (5): rampup upper bound constraint
                model.add_user_cut_constraints((x_var[g,t] <= uc_start_ub[gen_type]*y_var[g,t] + (uc_gen_ub[gen_type] - uc_start_ub[gen_type])*\
                    (y_var[g,t+1] - u_var[g,t+1] - u_var[g,t]) for g in range(start_index+1, end_index+1) for t in range(2, uc_opt_period)), \
                    names = 'three_period_conv_5')
                    
                #constraint (6): rampup upper bound constraint
                model.add_user_cut_constraints((x_var[g,t] <= uc_gen_ub[gen_type]*y_var[g,t] - (uc_gen_ub[gen_type] - uc_start_ub[gen_type])*\
                    u_var[g,t] - (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - uc_rampup_ub[gen_type])*u_var[g,t-1] \
                    for g in range(start_index+1, end_index+1) for t in range(3, uc_opt_period+1)), names = 'three_period_conv_6')
                    
                #constraint (7): two variables rampup upper bound constraint
                model.add_user_cut_constraints((x_var[g,t] - x_var[g,t-1] <= uc_start_ub[gen_type]*y_var[g, t] - uc_gen_lb[gen_type]*y_var[g,t-1] +\
                    (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*(y_var[g,t+1] - u_var[g,t+1] - u_var[g,t]) \
                    for g in range(start_index+1, end_index+1) for t in range(2, uc_opt_period)), names = 'three_period_conv_7')
                
                #constraint (10): two variables rampup upper bound constraint
                model.add_user_cut_constraints((x_var[g,t] - x_var[g,t+1] <= uc_start_ub[gen_type]*y_var[g,t] - uc_gen_lb[gen_type]*y_var[g,t+1] +\
                    (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*(y_var[g,t+1] - u_var[g,t+1] - u_var[g,t]) \
                    for g in range(start_index+1, end_index+1) for t in range(2,uc_opt_period)), names = 'three_period_conv_10')
                    
                #constraint (11): two variables rampup upper bound constraint
                model.add_user_cut_constraints((x_var[g,t] - x_var[g,t-2] <= (uc_gen_lb[gen_type]+2*uc_rampup_ub[gen_type])*y_var[g,t] -\
                    uc_gen_lb[gen_type]*y_var[g,t-2] - (uc_gen_lb[gen_type] + 2*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*u_var[g,t] -\
                    (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*u_var[g,t-1] for g in range(start_index+1, end_index+1) \
                    for t in range(3, uc_opt_period+1)), names = 'three_period_conv_11')
                    
                #constraint (12): two variables rampup upper bound constraint
                model.add_user_cut_constraints((x_var[g,t] - x_var[g,t+2] <= uc_start_ub[gen_type]*y_var[g,t] - uc_gen_lb[gen_type]*y_var[g,t+2] +\
                    uc_rampup_ub[gen_type]*(y_var[g,t+1] - u_var[g,t+1]) + (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*\
                    (y_var[g,t+2] - u_var[g,t+2] - u_var[g,t+1]) for g in range(start_index+1, end_index+1) for t in range(1, uc_opt_period-1)), \
                    names = 'three_period_conv_12')
                    
                #constraint (13): three variables rampup upper bound constraint
                model.add_user_cut_constraints((x_var[g,t] - x_var[g,t+1] + x_var[g,t+2] <= uc_start_ub[gen_type]*y_var[g,t] - (uc_start_ub[gen_type]-\
                    uc_rampup_ub[gen_type])*y_var[g,t+1] + uc_start_ub[gen_type]*y_var[g,t+2] + (uc_gen_ub[gen_type] - uc_start_ub[gen_type])*\
                    (y_var[g,t+2] - u_var[g,t+2] - u_var[g,t+1]) for g in range(start_index+1, end_index+1) for t in range(1, uc_opt_period-1)), \
                    names = 'three_period_conv_13')
            '''        
            '''add the 3period convex hull cuts with L=1=1 for L >= 1 of 2bin model to 2bin model as user cuts'''
            '''
            #constraint (24d):rampup upper bound constraint
            model.add_user_cut_constraints((x_var[g, t] <= uc_start_ub[gen_type]*y_var[g, t] + uc_rampup_ub[gen_type]*(y_var[g, t+1] - \
                u_var[g, t+1]) + (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - uc_rampup_ub[gen_type])*(y_var[g, t+2] - u_var[g, \
                t+2]) for g in range(start_index+1, end_index+1) for t in range(1, uc_opt_period-1)), names = 'three_period_conv_24d')
                                                                                                                               
            #constraint (24h): rampup upper bound constraint
            model.add_user_cut_constraints((x_var[g,t] <= (uc_start_ub[gen_type]+uc_rampup_ub[gen_type])*y_var[g,t] - \
                uc_rampup_ub[gen_type]*u_var[g,t] + (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - uc_rampup_ub[gen_type])*\
                (y_var[g, t-1] - u_var[g, t-1]) for g in range(start_index+1, end_index+1) for t in range(3, uc_opt_period+1)),\
                names = 'three_period_conv_24h')
            
            #constraint (24i): rampup upper bound constraint
            model.add_user_cut_constraints((x_var[g, t] - x_var[g, t-1] <= uc_start_ub[gen_type]*y_var[g,t] - \
                uc_gen_lb[gen_type]*y_var[g,t-1] + (uc_gen_lb[gen_type]+uc_rampup_ub[gen_type]-uc_start_ub[gen_type])*\
                (y_var[g,t+1] - u_var[g, t+1]) for g in range(start_index+1, end_index+1) for t in range(2, uc_opt_period)),\
                names =  'three_period_conv_24i')
                
            #constraint (24o): rampup upper bound constriant
            model.add_user_cut_constraints((x_var[g,t] - x_var[g, t-2] <= (uc_gen_lb[gen_type] + 2*uc_rampup_ub[gen_type])*y_var[g,t] -\
                uc_gen_lb[gen_type]*y_var[g,t-2] - (uc_gen_lb[gen_type] + 2*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*\
                u_var[g,t] for g in range(start_index+1, end_index+1) for t in range(3, uc_opt_period+1)), \
                names = 'three_period_conv_24o')
                
            #constraint (24p): rampup upper bound constriant
            model.add_user_cut_constraints((x_var[g, t] - x_var[g, t-2] <= (uc_start_ub[gen_type] + uc_rampup_ub[gen_type])*y_var[g, t] -\
                uc_rampup_ub[gen_type]*u_var[g,t] - uc_gen_lb[gen_type]*y_var[g, t-2] + (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type]-\
                uc_start_ub[gen_type])*(y_var[g,t-1] - u_var[g,t-1]) for g in range(start_index+1, end_index+1) for t in range(3, \
                uc_opt_period+1)), names = 'three_period_conv_24p')
            
            #constraint (24q): rampup upper bound constraint
            model.add_user_cut_constraints((x_var[g, t] - x_var[g, t+2] <= uc_start_ub[gen_type]*y_var[g,t] - uc_gen_lb[gen_type]*\
                y_var[g, t+2] + (uc_gen_lb[gen_type] + 2*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*(y_var[g, t+1] -\
                u_var[g, t+1]) for g in range(start_index+1, end_index+1) for t in range(1, uc_opt_period-1)), \
                names = 'three_period_conv_24q')
                
            #constraint (24r): rampup upper bound constraint
            model.add_user_cut_constraints((x_var[g,t] - x_var[g,t+2] <= uc_start_ub[gen_type]*y_var[g,t] - uc_gen_lb[gen_type]*y_var[g,t+2] +\
                uc_rampup_ub[gen_type]*(y_var[g, t+1] - u_var[g, t+1]) + (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] -\
                uc_start_ub[gen_type])*(y_var[g,t+2] - u_var[g,t+2]) for g in range(start_index+1, end_index+1) \
                for t in range(1, uc_opt_period-1)), names = 'three_period_conv_24r')
            '''    
            '''add multi-period 2-binary usecuts'''
            '''
            #2bin constraint (28): single var ramp-up upper bound cut
            ub = min(uc_min_on[gen_type], math.floor((uc_gen_ub[gen_type] - uc_start_ub[gen_type])/uc_rampup_ub[gen_type])+2)
            for k in range(1, ub+1):
                model.add_user_cut_constraints((x_var[g, t] <= uc_start_ub[gen_type]*y_var[g, t] + (uc_gen_ub[gen_type] - \
                     uc_start_ub[gen_type])*(y_var[g,t+1] - u_var[g, t+1]) - sum((uc_gen_ub[gen_type] -\
                     uc_start_ub[gen_type] - (s-1)*uc_rampup_ub[gen_type])*u_var[g, t-s+1] for s in range(1, k)) for g in \
                     range(start_index+1, end_index+1) for t in range(k, uc_opt_period)), names = '2bin_usercuts_28')
            
            #2bin constraint (29): single variable rampup upper bound constraint
            #here we only consider the case that S includes all elements between t-m+1 and t, that is S = {t-m+1, t-m+2,...,t}
            for t in range(uc_min_on[gen_type]+1, uc_opt_period+1):
                #define the upper bound of m
                m_ub = min(t-uc_min_on[gen_type]-1, max(math.floor((uc_gen_ub[gen_type]-uc_start_ub[gen_type])/uc_rampup_ub[gen_type])-\
                    uc_min_on[gen_type]+1, 0))
                for m in range(0, m_ub+1):
                    if t-m+1 <= t:
                        #we should ensure that S is not empty set
                        #calculate all parameter values
                        beta_t = uc_gen_ub[gen_type] - uc_start_ub[gen_type] - (m+uc_min_on[gen_type]-1)*uc_rampup_ub[gen_type]
                        # we add cut (29) to each m in this range
                        n_lb = min(1, uc_min_on[gen_type]-1, uc_opt_period-t)
                        n_ub = min(uc_min_on[gen_type]-1, uc_opt_period-t)
                        for n in range(n_lb, n_ub+1):
                            #use n_round to denote [n-1]^+
                            n_round = max(n-1, 0)
                            if n <= uc_opt_period-t-1 and n < (uc_min_on[gen_type]-1)/2:
                                #if n <= T-t-1 and n < (L-1)/2, we do not add the cut to the model because it is not valid and facet-defining
                                pass
                            else:
                                model.add_user_cut_constraints((x_var[g,t] <= uc_start_ub[gen_type]*y_var[g,t] +\
                                    uc_rampup_ub[gen_type]*sum(y_var[g, i] - model.sum_vars(u_var[g,i-j] for j in range(0, uc_min_on[gen_type])) \
                                    for i in range(t-m+1, t+1)) + uc_rampup_ub[gen_type]*sum(y_var[g,t+k] - model.sum_vars(u_var[g, t+k-j] \
                                    for j in range(0, uc_min_on[gen_type])) for k in range(1, n_round+1)) + (uc_min_on[gen_type]-1-n_round)*\
                                    uc_rampup_ub[gen_type]*(y_var[g, t+n] - model.sum_vars(u_var[g, t+n-j] for j in range(0, uc_min_on[gen_type])))+\
                                    beta_t*(y_var[g, t-m] - model.sum_vars(u_var[g, t-m-j] for j in range(0, uc_min_on[gen_type]))) +\
                                    uc_rampup_ub[gen_type]*sum(k*u_var[g, t-k] for k in range(1, t+uc_min_on[gen_type]-uc_opt_period)) +\
                                    uc_rampup_ub[gen_type]*sum(min(uc_min_on[gen_type]-1-k, k)*u_var[g, t-k] for k in range(max(t+uc_min_on[gen_type]-\
                                    uc_opt_period,0), uc_min_on[gen_type])) for g in range(start_index+1, end_index+1)), names = '2bin_usercut_29')
            
            #2bin constraint (30): single variable rampup uuper bound constraint
            #we only consider the case that S includes all elements between t_hat+1 and t+m, that is S = {t_hat+1,...,t+m}
            for t in range(1, uc_opt_period):
                if min(t-2, uc_min_on[gen_type]-2) >= uc_min_on[gen_type]/2:
                    t_hat = t + min(t-2, uc_min_on[gen_type]-2)
                else:
                    t_hat = max(t+1, uc_min_on[gen_type]+1)
                #define the lower and upper bound for m
                m_lb = max(t_hat-t-1, 0)
                m_ub = min(uc_opt_period-t-1, math.floor((uc_gen_ub[gen_type] - uc_start_ub[gen_type])/uc_rampup_ub[gen_type]))
                for m in range(m_lb, m_ub+1):
                    if t_hat <= t+m:
                        #we should ensure S is not empty, and it could be empty if t_hat = t+m
                        model.add_user_cut_constraints((x_var[g, t] <= uc_start_ub[gen_type]*y_var[g,t] + uc_rampup_ub[gen_type]*sum(y_var[g,i] -\
                            model.sum_vars(u_var[g, i-j] for j in range(0, min(uc_min_on[gen_type]-1, i-2)+1)) for i in range(t+1, t_hat)) +\
                            uc_rampup_ub[gen_type]*sum(y_var[g,i] - model.sum_vars(u_var[g, i-j] for j in range(0, uc_min_on[gen_type])) \
                            for i in range(t_hat, t+m+1)) + (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - m*uc_rampup_ub[gen_type])*(y_var[g, t+m+1] -\
                            model.sum_vars(u_var[g, t+m+1-j] for j in range(0, uc_min_on[gen_type]))) + uc_rampup_ub[gen_type]*sum(k*u_var[g,t-k] \
                            for k in range(1, t+uc_min_on[gen_type]-uc_opt_period)) + uc_rampup_ub[gen_type]*sum(min(uc_min_on[gen_type]-1-k, k)*\
                            u_var[g, t-k] for k in range(max(t+uc_min_on[gen_type]-uc_opt_period, 0), min(uc_min_on[gen_type]-1, t-2)+1)) \
                            for g in range(start_index+1, end_index+1)), names = '2bin_usercut_30')
            
            #2bin constraint (31) & (32): two variables rampup upper bound constraint
            #For cut (32), we only consider the that case S includes all elements between t-m+L and t
            for t in range(2, uc_opt_period+1):
                #define the upper bound of m
                m_ub = min(t-1, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampup_ub[gen_type]))
                for m in range(1, m_ub+1):
                    #define the lower and upper bound for n
                    n_lb = min(1, uc_opt_period-t)
                    n_ub = min(m, uc_min_on[gen_type], uc_opt_period-t)
                    for n in range(n_lb, n_ub+1):
                        #valid and facet-defining condition (i): if min{m-1, L-2} >= L/2 then n >= min{m-1, L-2, T-t} 
                        if min(m-1, uc_min_on[gen_type]-2) >= uc_min_on[gen_type]/2 and n < min(m-1, uc_min_on[gen_type]-2, uc_opt_period-t):
                            pass
                        else:
                            if m <= uc_min_on[gen_type]-1:
                                #if S is an empty set, then if n <= L-1-m, we require n >= min(m, L, T-t)
                                if n <= uc_min_on[gen_type]-1-m and n < min(m, uc_min_on[gen_type], uc_opt_period-t):
                                    pass
                                else:
                                    #constraint (31)
                                    model.add_user_cut_constraints((x_var[g, t] - x_var[g, t-m] <= uc_start_ub[gen_type]*y_var[g, t] -\
                                        uc_gen_lb[gen_type]*y_var[g,t-m] + uc_rampup_ub[gen_type]*sum(y_var[g, t+k] - model.sum_vars(u_var[g,t+k-j] \
                                        for j in range(0, min(uc_min_on[gen_type], k+m))) for k in range(1, max(n-1, 0)+1)) + (uc_gen_lb[gen_type] +\
                                        (m-max(n-1, 0))*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*(y_var[g,t+n] - model.sum_vars(u_var[g, t+n-j] \
                                        for j in range(0, min(uc_min_on[gen_type], n+m)))) + uc_rampup_ub[gen_type]*sum(k*u_var[g, t-k] \
                                        for k in range(1, min(t+uc_min_on[gen_type]-uc_opt_period, m))) + uc_rampup_ub[gen_type]*\
                                        sum(min(uc_min_on[gen_type]-1-k, k)*u_var[g, t-k] for k in range(max(t+uc_min_on[gen_type]-uc_opt_period, 0),\
                                        min(uc_min_on[gen_type], m))) for g in range(start_index+1, end_index+1)), names = '2bin_usercut_31')
                            else:
                                pass
            '''
            '''
                            #if S is not empty, when min{m-1, L-2} >= L/2, n >= min{m-1, L-2, T-t}
                            if min(m-1, uc_min_on[gen_type]-2) >= uc_min_on[gen_type]/2 and n < min(m-1, uc_min_on[gen_type]-2, uc_opt_period-t):
                                pass
                            else:
                                #since here q is not given and d_i is unknown, so we omit cut (32) here!
                                #constraint (32)
                                model.add_user_cut_constraints((x_var[g,t] - x_var[g,t-m] <= uc_start_ub[gen_type]*y_var[g,t] -\
                                    uc_gen_lb[gen_type]*y_var[g,t-m] + uc_rampup_ub[gen_type]*sum(y_var[g,i] - model.sum_vars(u_var[g,i-j] \
                                    for j in range(0, uc_min_on[gen_type])) for i in range(t-m+uc_min_on[gen_type]+1, t+1)) + uc_rampup_ub[gen_type]*\
                                    sum(y_var[g,t+k] - model.sum_vars(u_var[g, t+k-j] for j in range(0, min(uc_min_on[gen_type], k+m))) \
                                    for k in range(1, max(n-1, 0)+1)) + (uc_min_on[gen_type]-1-max(n-1, 0))*uc_rampup_ub[gen_type]*(y_var[g,t+n] -\
                                    model.sum_vars(u_var[g,t+n-j] for j in range(0, min(uc_min_on[gen_type], n+m)))) + (uc_gen_lb[gen_type] +\
                                    uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*(y_var[g,q] - model.sum_vars(u_var[g,q-j] \
                                    for j in range(0, min(uc_min_on[gen_type], q-t+m)))) + uc_rampup_ub[gen_type]*sum(k*u_var[g, t-k] \
                                    for k in range(1, min(t+uc_min_on[gen_type]-uc_opt_period, m))) + uc_rampup_ub[gen_type]*\
                                    sum(min(uc_min_on[gen_type]-1-k, k)*u_var[g, t-k] for k in range(max(t+uc_min_on[gen_type]-uc_opt_period, 0),\
                                    min(uc_min_on[gen_type], m))) for g in range(start_index+1, end_index+1)), names = '2bin_usercut_32')
            '''
            '''
            #2bin constraint (33): two vairables rampup upper bound
            for t in range(1, uc_opt_period):
                #define parameters for calculating
                if min(t-2, uc_min_on[gen_type]-2) >= uc_min_on[gen_type]/2:
                    t_hat = t + min(t-2, uc_min_on[gen_type]-2)
                else:
                    t_hat = max(t+1, uc_min_on[gen_type]+1)
                #define the upper bound for m
                m_ub = min(uc_opt_period-t, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampup_ub[gen_type]))
                for m in range(1, m_ub+1):
                    t_tilde = min(t_hat, t+m)
                    if t+1 <= t_tilde-1 and t_tilde+1 <= t+m:
                        #if S_0, S_1, S is not empty
                        q = t+m
                        if m <= uc_min_on[gen_type]-1 and m < math.floor((uc_min_on[gen_type]+1)/2):
                            #since q = t+m, when m <= L-1, m >= math.floor((L+1)/2)
                            pass
                        else:
                            #add the cut (33)
                            model.add_user_cut_constraints((x_var[g,t] - x_var[g, t+m] <= uc_start_ub[gen_type]*y_var[g,t] - uc_gen_lb[gen_type]*\
                                y_var[g, t+m] + uc_rampup_ub[gen_type]*sum(y_var[g,i] - model.sum_vars(u_var[g,i-j] \
                                for j in range(0, min(uc_min_on[gen_type], i-1))) for i in range(t+1, t_tilde)) + uc_rampup_ub[gen_type]*\
                                sum(y_var[g,i] - model.sum_vars(u_var[g,i-j] for j in range(0, uc_min_on[gen_type])) for i in range(t_hat, t+m)) +\
                                (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*(y_var[g,q] - model.sum_vars(u_var[g, q-j] \
                                for j in range(0, min(uc_min_on[gen_type], q-1)))) + uc_rampup_ub[gen_type]*sum(k*u_var[g, t-k] \
                                for k in range(1, t+uc_min_on[gen_type]-uc_opt_period)) + uc_rampup_ub[gen_type]*sum(min(uc_min_on[gen_type]-1-k, k)*\
                                u_var[g,t-k] for k in range(max(t+uc_min_on[gen_type]-uc_opt_period, 0), min(uc_min_on[gen_type], t-1))) \
                                for g in range(start_index+1, end_index+1)), names = '2bin_usercut_33')   
                                                                   
            #2bin constraint (34) for 2 bin model: double var ramp-up upper bound cut
            if uc_min_on[gen_type] >= 2:
                #the cut is only valid for L >= 2
                t_lb = max(uc_min_on[gen_type]+1, 3)
                model.add_user_cut_constraints((x_var[g, t-2] - x_var[g, t-1] + x_var[g, t] <= uc_start_ub[gen_type]*y_var[g, t-2] - \
                    (uc_start_ub[gen_type] - uc_rampup_ub[gen_type])*y_var[g, t-1] + uc_start_ub[gen_type]*y_var[g, t] + \
                    (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*(y_var[g, t+1] - u_var[g, t+1] - \
                    y_var[g, t]) + (uc_gen_ub[gen_type] - uc_start_ub[gen_type])*(y_var[g,t] - u_var[g,t] - u_var[g, t-1]) - \
                    sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - s*uc_rampup_ub[gen_type])*u_var[g, t-s-2] for s in \
                    range(0, uc_min_on[gen_type]-2)) for g in range(start_index+1, end_index+1) for t in range(t_lb, uc_opt_period)), \
                    names = '2bin_usercuts_34')

            #2bin constaint (35): three variables rampup upper bound
            if uc_min_on[gen_type] >= 2:
                #the cut is valid for L >= 2
                for t in range(uc_min_on[gen_type]+1, uc_opt_period+1):
                    #define the lower and upper bound for m
                    m_lb = max(3-uc_min_on[gen_type], 0)
                    m_ub = min(max(t-uc_min_on[gen_type]-1, 0), max(math.floor((uc_gen_ub[gen_type] -\
                           uc_start_ub[gen_type])/uc_rampup_ub[gen_type])-uc_min_on[gen_type]+3, 0))
                    for m in range(m_lb, m_ub+1):
                        beta_t = uc_gen_ub[gen_type] - uc_start_ub[gen_type] - uc_rampup_ub[gen_type]*max(m+uc_min_on[gen_type]-3, 0)
                        if uc_min_on[gen_type] == 2:
                            #if L = 2, we add the cut, which takes different form with L >= 3
                            model.add_user_cut_constraints((x_var[g,t-2] - x_var[g,t-1] + x_var[g,t] <= uc_start_ub[gen_type]*y_var[g,t-2] -\
                                (uc_start_ub[gen_type] - uc_rampup_ub[gen_type])*y_var[g,t-1] + uc_start_ub[gen_type]*y_var[g,t] + uc_rampup_ub[gen_type]*\
                                sum(y_var[g,i] - model.sum_vars(u_var[g,i-j] for j in range(0,uc_min_on[gen_type])) for i in range(t-m+1, t-1))-\
                                max(uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - uc_start_ub[gen_type], 0)*u_var[g,t-2] +\
                                (max(m+uc_min_on[gen_type]-3, 0)-m+2)*uc_rampup_ub[gen_type]*(y_var[g,t] - model.sum_vars(u_var[g,t-k] \
                                for k in range(0, uc_min_on[gen_type]))) + sum((k-2)*uc_rampup_ub[gen_type]*u_var[g,t-k] \
                                for k in range(3, uc_min_on[gen_type])) + beta_t*(y_var[g,t-m] - model.sum_vars(u_var[g,t-m-j] \
                                for j in range(0, uc_min_on[gen_type]))) for g in range(start_index+1, end_index+1)), names = '2bin_usercut_35_I')
                        else:
                            #if L >= 3, cut takes another form as below
                            model.add_user_cut_constraints((x_var[g,t-2] - x_var[g,t-1] + x_var[g,t] <= uc_start_ub[gen_type]*y_var[g,t-2] -\
                                (uc_start_ub[gen_type] - uc_rampup_ub[gen_type])*y_var[g,t-1] + uc_start_ub[gen_type]*y_var[g,t] + uc_rampup_ub[gen_type]*\
                                sum(y_var[g,i] - model.sum_vars(u_var[g,i-j] for j in range(0,uc_min_on[gen_type])) for i in range(t-m+1, t)) +\
                                (max(m+uc_min_on[gen_type]-3, 0)-m+1)*uc_rampup_ub[gen_type]*(y_var[g,t] - model.sum_vars(u_var[g,t-k] \
                                for k in range(0, uc_min_on[gen_type]))) + sum((k-2)*uc_rampup_ub[gen_type]*u_var[g,t-k] \
                                for k in range(3, uc_min_on[gen_type])) + beta_t*(y_var[g,t-m] - model.sum_vars(u_var[g,t-m-j] \
                                for j in range(0, uc_min_on[gen_type]))) for g in range(start_index+1, end_index+1)), names = '2bin_usercut_35_II')
            
            #2bin constraint (36): three variables rampup upper bound
            for t in range(3, uc_opt_period):
                #define all parameters
                t_hat = max(t+1, uc_min_on[gen_type]+1)
                #define the lower and upper bound for m
                m_lb = max(t_hat-t-1, 0)
                m_ub = min(uc_opt_period-t-1, math.floor((uc_gen_ub[gen_type] - uc_start_ub[gen_type])/uc_rampup_ub[gen_type]))
                for m in range(m_lb, m_ub+1):
                    if t+1 <= t_hat-1 and t_hat <= t+m:
                        #S_1 can not be empty and S can be empty but t_hat <= t+m
                        #add the cut (36)
                        model.add_user_cut_constraints((x_var[g,t-2] - x_var[g,t-1] + x_var[g,t] <= uc_start_ub[gen_type]*y_var[g,t-2] -\
                            (uc_start_ub[gen_type] - uc_rampup_ub[gen_type])*y_var[g,t-1] + uc_start_ub[gen_type]*y_var[g,t] + uc_rampup_ub[gen_type]*\
                            sum(y_var[g,i] - model.sum_vars(u_var[g,i-j] for j in range(0, min(uc_min_on[gen_type], i-1))) for i in range(t+1, t_hat))+\
                            uc_rampup_ub[gen_type]*sum((k-2)*u_var[g,t-k] for k in range(3, min(t-1, uc_min_on[gen_type]))) + uc_rampup_ub[gen_type]*\
                            sum(y_var[g,i] - model.sum_vars(u_var[g,i-j] for j in range(0, uc_min_on[gen_type])) for i in range(t_hat, t+m+1)) +\
                            (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - m*uc_rampup_ub[gen_type])*(y_var[g,t+m+1] - model.sum_vars(u_var[g,t+m+1-j] \
                            for j in range(0, uc_min_on[gen_type]))) for g in range(start_index+1, end_index+1)), names = '2bin_usercut_36')
            '''
            start_index = end_index
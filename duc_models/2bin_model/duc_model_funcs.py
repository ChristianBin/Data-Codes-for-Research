# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 09:36:14 2022

@author: tianb
"""

'''This code is for solving a 2-bin UC model'''
from docplex.mp.callbacks.cb_mixin import *
from docplex.mp.model import Model
from duc_callback import *
from duc_params import *
from duc_export_results import *
from duc_cutpool import *

def set_variables(duc_model, instance_index):
    '''adding decision variables for the model'''
    
    #here because in main function the instance_index starts from 1 to 20, we use instance_index - 1. 
    #setting generation amount decision varibale x[g][t] for each generator g in [1, G] for each period t in [1, T]
    x_var = {(g,t): duc_model.continuous_var(name = 'x_{0}_{1}'.format(g, t), lb = 0) \
             for g in range(1, uc_total_gen_num[instance_index-1] + 1) for t in range(1, uc_opt_period + 1)}
        
    #setting on/off status decision variable y[g][t] for each generator g in [1, G] for each period t in [1, T]
    y_var = {(g,t): duc_model.binary_var(name = 'y_{0}_{1}'.format(g,t)) \
             for g in range(1, uc_total_gen_num[instance_index-1] + 1) for t in range(1, uc_opt_period + 1)}
    
    #setting start-up/not status decision variable u[g,t] for each generator in [1, G] and each period t in [2,T]
    u_var = {(g,t): duc_model.binary_var(name = 'u_{0}_{1}'.format(g,t)) \
             for g in range(1, uc_total_gen_num[instance_index-1] + 1) for t in range(2, uc_opt_period +1)}
        
    #setting generation cost decision variable for each generator g in [1, G] for each period t in [1, T]
    fx_var = {(g,t): duc_model.continuous_var(name = 'fx_{0}_{1}'.format(g,t), lb = 0) \
              for g in range(1, uc_total_gen_num[instance_index-1] + 1) for t in range(1, uc_opt_period + 1)}
    
    return x_var, y_var, u_var, fx_var

def add_polytope_constraints(duc_model, x_var, y_var, u_var, fx_var, instance_index):
    '''adding polytope constraints to the model'''
        
    #we use start_index and end_index to denote the index of current type of generator in the generator
    #array [1, G]
    start_index = 0
    end_index = 0
    #add all constraints for one type of generators at each iteration, totally 8 types
    for gen_type in range(0, 8):
        #if the number of the current type of generator is not 0, we add the minimum on constriant
        if uc_gen_num[instance_index-1, gen_type] > 0:
            #update the end_index of this type of generator in [1, G]
            end_index += uc_gen_num[instance_index-1, gen_type]
            
            #constraint (1a) and (1b): minimum on and off period constraints for each generator g in [1, G]
            for g in range(start_index+1, end_index+1):
                for t in range(uc_min_on[gen_type]+1, uc_opt_period+1):
                    duc_model.add_constraint((duc_model.sum_vars(u_var[g, i] for i in range(t-uc_min_on[gen_type]+1, t+1)) <= \
                        y_var[g, t]), ctname = 'cstr_min_on_{0}_{1}'.format(g,t))
                        
                for t in range(uc_min_off[gen_type]+1, uc_opt_period+1):
                    duc_model.add_constraint((duc_model.sum_vars(u_var[g,i] for i in range(t-uc_min_off[gen_type]+1, t+1)) <= \
                        1 - y_var[g, t-uc_min_off[gen_type]]), ctname = 'cstr_min_off_{0}_{1}'.format(g, t))
            
            #constraint (1c): relationship between u_var and y_var
            for g in range(start_index+1, end_index+1):
                for t in range(2, uc_opt_period+1):
                    duc_model.add_constraint((y_var[g, t] - y_var[g, t-1] - u_var[g,t] <= 0), ctname = 'logic_y&u_{0}_{1}'.format(g,t))
            
            #constraint (1d): generation lower bound
            duc_model.add_constraints((-x_var[g,t] + uc_gen_lb[gen_type]*y_var[g,t] <= 0 \
                                       for t in range(1, uc_opt_period + 1) for g in range(start_index + 1, \
                                       end_index + 1)), names = 'cstr_gen_lb')
            #constraint (1e): generation upper bound
            duc_model.add_constraints((x_var[g,t] - uc_gen_ub[gen_type]*y_var[g,t] <= 0 \
                                       for t in range(1, uc_opt_period +1) for g in range(start_index + 1, \
                                       end_index + 1)), names = 'cstr_gen_ub')
            
            #constraint(1f): ramp up limitation for successive period for each type of generator
            duc_model.add_constraints((x_var[g,t] - x_var[g,t-1] <= uc_rampup_ub[gen_type]*y_var[g,t-1] + \
                                       uc_start_ub[gen_type]*(1 - y_var[g,t-1]) for t in \
                                       range(2, uc_opt_period + 1) for g in \
                                       range(start_index + 1, end_index + 1)), names = 'cstr_rampup')
            #constraint (1g): ramp down limitation for successive period for each type of generator
            duc_model.add_constraints((x_var[g,t-1] - x_var[g,t] <= uc_rampdown_ub[gen_type]*y_var[g,t] + \
                                      uc_start_ub[gen_type]*(1 - y_var[g,t]) for t in \
                                      range(2, uc_opt_period + 1) for g in \
                                      range(start_index + 1, end_index + 1)), names ='cstr_rampdown')
            
            #constraint (1h): generation cost constraint fx[g][t] in each period [1, T] for each generator
            #[1, G], we  set it as a decision varible
            #the length of one part if we divide [uc_gen_lb[i], uc_gen_ub[i]] into uc_gen_cost_pwl_num[i] 
            #parts
            unit_part_length = (uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_gen_cost_pwl_num[gen_type]
            for k in range(0, uc_gen_cost_pwl_num[gen_type] + 1):
                #we break x_var[g][t] into uc_gen_cost_pwl_num[i] parts and then calculate the slope of
                #f(x) = ax^2 + bx + c at that point to get the linear expression of its tangent
                #then we let fx_var[g][t] is greater than all the tangents at these points to approximate
                #f(x). E.g. at point x1, the tangent line expression is t(x) = (2ax1 + b)(x - x1) + 
                #ax1^2 + bx^1 + c
                #x_k, the piecewise point we choose
                piecewise_point = uc_gen_lb[gen_type] + k*unit_part_length
                duc_model.add_constraints((fx_var[g,t] >= (2*uc_gen_cost_func_a[gen_type]*piecewise_point + \
                    uc_gen_cost_func_b[gen_type])*x_var[g,t] - (uc_gen_cost_func_a[gen_type]*piecewise_point**2 \
                    - uc_gen_cost_func_c[gen_type])*y_var[g,t] for t in range(1, uc_opt_period + 1) for \
                    g in range(start_index + 1, end_index + 1)), names = 'fx_pwl_cost')
                                 
            '''add 1bin two-period convex hull cuts to 2bin model'''
            
            #constraint (2a): upper bound of y
            #duc_model.add_constraints((y_var[g, t] <= 1 for g in range(start_index+1, end_index+1) for t in range(1, \
                #uc_opt_period+1)), names = 'two_period_conv_a')
            #constraint (2c): start up upper bound constraint
            duc_model.add_constraints((x_var[g, t] <= uc_start_ub[gen_type]*y_var[g, t] + (uc_gen_ub[gen_type] - \
                uc_start_ub[gen_type])*y_var[g, t+1] for g in range(start_index+1, end_index+1) for t in range(1, \
                uc_opt_period)), names = 'two_period_conv_c')
            #constraint (2d): shut down upper bound constraint
            duc_model.add_constraints((x_var[g, t] <= (uc_gen_ub[gen_type] - uc_start_ub[gen_type])*y_var[g, t-1] + \
                uc_start_ub[gen_type]*y_var[g,t] for g in range(start_index+1, end_index+1) for t in range(2, \
                uc_opt_period+1)), names = 'two_period_conv_d')
            #constraint (2e): ramp-up upper bound constraint
            duc_model.add_constraints((x_var[g, t] - x_var[g, t-1] <= (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type])*y_var[g,t]- \
                uc_gen_lb[gen_type]*y_var[g, t-1] for g in range(start_index+1, end_index+1) for t in range(2, uc_opt_period + 1)), \
                names = 'two_period_conv_e')
            #constraint (2f): start-up upper bound constraint
            duc_model.add_constraints((x_var[g,t] - x_var[g, t-1] <= uc_start_ub[gen_type]*y_var[g, t] - (uc_start_ub[gen_type] - \
                uc_rampup_ub[gen_type])*y_var[g, t-1] for g in range(start_index+1, end_index+1) for t in range(2, uc_opt_period+1)), \
                names = 'two_period_conv_f')
            #constraint (2g): ramp-down upper bound constraint
            duc_model.add_constraints((x_var[g, t] - x_var[g, t+1] <= (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type])*y_var[g, t] - \
                uc_gen_lb[gen_type]*y_var[g, t+1] for g in range(start_index+1, end_index+1) for t in range(1, uc_opt_period)), \
                names = 'two_period_conv_g')
            #constraint (2h): shut-down upper bound constraint
            duc_model.add_constraints((x_var[g, t] - x_var[g, t+1] <= uc_start_ub[gen_type]*y_var[g, t] - (uc_start_ub[gen_type] - \
                uc_rampup_ub[gen_type])*y_var[g, t+1] for g in range(start_index+1, end_index+1) for t in range(1, uc_opt_period)), \
                names = 'two_period_conv_h')          
            
            
            '''add 2bin two-period convex hull cuts to 2bin model as constraints'''
            '''
            #since (2a)-(2c) have been included in the original polytope, here we only add (2d)-(2g) to the 2bin model
            #constraint (2d): start-up upper bound constraint
            duc_model.add_constraints((x_var[g,t] <= uc_start_ub[gen_type]*y_var[g, t] + (uc_gen_ub[gen_type] - \
                uc_start_ub[gen_type])*(y_var[g, t+1] - u_var[g, t+1]) for g in range(start_index+1, end_index+1) \
                for t in range(1, uc_opt_period)), names = 'two_period_conv_d')
            
            #constraint(2e): start-up upper bound constraint
            duc_model.add_constraints((x_var[g, t] <= uc_gen_ub[gen_type]*y_var[g,t] - (uc_gen_ub[gen_type] -\
                uc_start_ub[gen_type])*u_var[g,t] for g in range(start_index+1, end_index+1) for t in range(2, \
                uc_opt_period+1)), names = 'two_period_conv_e')
                                                                                                            
            #constraint (2f): ramp-up uppper bound constraint
            duc_model.add_constraints((x_var[g, t+1] - x_var[g, t] <= (uc_gen_lb[gen_type]+uc_rampup_ub[gen_type])*y_var[g, t+1]- \
                uc_gen_lb[gen_type]*y_var[g,t] - (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*u_var[g,t+1] \
                for g in range(start_index+1, end_index+1) for t in range(1, uc_opt_period)), names = 'two_period_conv_f')
                
            #constriant (2g): ramp-down upper bound constraint
            duc_model.add_constraints((x_var[g,t] - x_var[g, t+1] <= uc_start_ub[gen_type]*y_var[g,t] - (uc_start_ub[gen_type] -\
                uc_rampup_ub[gen_type])*y_var[g,t+1] - (uc_gen_lb[gen_type] + uc_rampup_ub[gen_type] - \
                uc_start_ub[gen_type])*u_var[g,t+1] for g in range(start_index+1, end_index+1) for t in range(1, uc_opt_period)), \
                names = 'two_period_conv_g')
            '''
                          
            #update the start_index for next operation
            start_index = end_index
            
    #constriant (1i): total generation amount for all generators at each period should be equal to
    #the load at that time period
    duc_model.add_constraints((sum(x_var[g,t] for g in range(1, uc_total_gen_num[instance_index-1] + 1)) == \
                uc_load[instance_index-1,t-1] for t in range(1, uc_opt_period+1)), names = 'cstr_load')
        
    #constraint (1j): spinning reserve constraint
    duc_model.add_constraints((sum(uc_gen_ub[i]*y_var[g,t] for i in range(0, 8) for g in \
                    range(uc_gen_accumu_num[instance_index-1, i]+1, uc_gen_accumu_num[instance_index-1, i+1]+1)) \
                    >= (1 + uc_spin_reserve[t-1])*uc_load[instance_index-1, t-1] for t in \
                    range(1, uc_opt_period+1)), names = 'cstr_spin_reserve')
        
def set_obj_function(duc_model, x_var, y_var, u_var, fx_var, instance_index):
    '''setting the objective function for the model'''
    #summing up the total start-up  and shut-down cost
    total_start_cost = sum(u_var[g, t]*uc_unit_startup_cost[gen_type] for gen_type in range(0,8) for g in \
                       range(uc_gen_accumu_num[instance_index-1, gen_type]+1, uc_gen_accumu_num[instance_index-1, gen_type+1]+1) \
                       for t in range(2, uc_opt_period+1))
    total_shut_cost = sum((y_var[g, t-1] - y_var[g, t] + u_var[g, t])*uc_unit_shutdown_cost[gen_type] for gen_type in range(0,8) \
                      for g in range(uc_gen_accumu_num[instance_index-1, gen_type]+1, uc_gen_accumu_num[instance_index-1, \
                      gen_type+1]+1) for t in range(2, uc_opt_period+1))
    total_start_shut_cost = total_start_cost + total_shut_cost
    
    total_generation_cost = duc_model.sum_vars(fx_var[g,t] for g in \
                            range(1, uc_total_gen_num[instance_index-1] + 1) for t in range(1, uc_opt_period+1))
      
    #adding kpi
    duc_model.add_kpi(total_start_shut_cost, 'total start up and shut down cost')
    duc_model.add_kpi(total_generation_cost, 'total genration cost')
     
    #minimize the total cost    
    duc_model.minimize(total_start_shut_cost + total_generation_cost)
    
def build_model(usercutpool, usercutcallback, instance_index):
    '''creating an mp model for the problem'''
    
    #denote the model by duc_model
    m = Model(name = 'duc_model')
    #setting decision varibles
    x_var, y_var, u_var, fx_var = set_variables(m, instance_index)
    #adding polytope constraints
    add_polytope_constraints(m, x_var, y_var, u_var, fx_var, instance_index)
    
    '''
    #if do not use lazycut callback, add the cut as normal constraint or via Model.add_lazy_constraint()
    if not lazycutcallback:
        #m.add_constraints()
        pass
    '''
    
    # setting the objective function
    set_obj_function(m, x_var, y_var, u_var, fx_var, instance_index)
    
    '''
    #if use lazycut callback, add this callback to the model
    if lazycutcallback:
        #register a lazy constraint callback
        ('* add lazy constriant callback')
        lazyct_cb = m.register_callback(LazyCallback)
        
        #store lazy callback inside the callback, as DOcplex object
        lazyct_cb.add_lazy_constraints()
        print('add lazy constraints callback with {0}'.format(len(lazyct_cb.cts)))
        m.lazy_callback = lazyct_cb
    '''
    
    if usercutpool:
        '''if usercutpool is true, we add usercuts to usercut pool'''
        add_cutpool_cuts(m, x_var, y_var, u_var, instance_index)
    
    #if use usercut callback, add this callback to the model
    if usercutcallback:
        #register usercut callback
        user_cb = m.register_callback(UserCallback)
        user_cb.instance_index = instance_index
        user_cb.x_var = x_var
        user_cb.y_var = y_var
        #user_cb.u_var = u_var
    return m

def solve_model(usercutpool, usercutcallback, instance_index, num_of_run):
    '''solve the model we built'''
    
    #build a duc model
    duc_model = build_model(usercutpool, usercutcallback, instance_index)
    
    # Tweak some CPLEX parameters so that CPLEX has a harder time to
    # solve the model and our cut separators can actually kick in.
    #turn off the presolve to let our usercut callback kick in
    #duc_model.parameters.preprocessing.presolve = 0
    #change the search method to branch and cut
    duc_model.parameters.mip.strategy.search = 1
    #set relative mip gap tolerance
    #duc_model.parameters.mip.tolerances.mipgap = 1e-7
    #set absolute mip gap tolerance
    #duc_model.parameters.mip.tolerances.absmipgap = 1e-7
    #params = m.parameters
    #set the smphasis on numerical precision
    #duc_model.parameters.emphasis.numerical = 1
    #set the feasibility tolerance
    #duc_model.parameters.simplex.tolerances.feasibility = 1e-9
    #setting the thread to be 1
    #duc_model.parameters.threads = 1
    #duc_model.parameters.mip.strategy.heuristicfreq = -1
    #duc_model.parameters.mip.cuts.mircut = -1
    #duc_model.parameters.mip.cuts.implied = -1
    #duc_model.parameters.mip.cuts.gomory = -1
    #duc_model.parameters.mip.cuts.flowcovers = -1
    #duc_model.parameters.mip.cuts.pathcut = -1
    #duc_model.parameters.mip.cuts.liftproj = -1
    #duc_model.parameters.mip.cuts.zerohalfcut = -1
    #duc_model.parameters.mip.cuts.cliques = -1
    #duc_model.parameters.mip.cuts.covers = -1
    
    #set time limit on solving
    duc_model.set_time_limit(time_limit)
    
    #set the filename for different optimization strategy
    log_filename = set_file_name(usercutpool, usercutcallback) + '_logfile' 
    model_filename = set_file_name(usercutpool, usercutcallback) + '_model'
    log_filename = filepath + '/instance{0}/'.format(instance_index) + log_filename + \
        '_instance{0}_{1}.txt'.format(instance_index, num_of_run)
    model_filename = model_filename + '_instance{0}_{1}'.format(instance_index, num_of_run)
    #export the model file as .lp file
    duc_model.export_as_lp(path = filepath + '/instance{0}/'.format(instance_index), basename = model_filename)
    
    #solve the model and export the log file
    duc_model_solve = duc_model.solve(log_output = log_filename)
    
    #assert that it is solved successfully otherwise raise exception
    assert duc_model_solve
    duc_model.report()
    export_sol_result(duc_model, usercutpool, usercutcallback, instance_index, num_of_run)
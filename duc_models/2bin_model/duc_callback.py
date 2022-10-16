# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 21:39:01 2022

@author: tianb
"""
'''usercut callback class'''
import math
from cplex.callbacks import UserCutCallback
from cplex.callbacks import LazyConstraintCallback
from docplex.mp.model import Model
from docplex.mp.callbacks.cb_mixin import *
from duc_params import *
'''
class CutpoolCallback(ConstraintCallbackMixin, UserCutCallback):
    #usercut callback class by doing separation through docplex
    
    def __init__(self, env):
       UserCutCallback.__init__(self, env)
       ConstraintCallbackMixin.__init__(self)
       self.eps = 1e-1
       #fraction of constraints to pop out in self.cts
       #self.frac_pop = 1
       #count the number of cuts added to the model during solving
       self.nb_cuts = 0
            
    def add_cut_constraint(self, ct):
        self.register_constraint(ct)
            
    @print_called('--> usercut callback called: #{0}')
    def __call__(self):
        #fetch variable solutioin values at this point
        sol = self.make_complete_solution()
        #fetch those constraints which are not satisfied
        #checking unsatisfied cut and add it to the subproblem
            
        unsats = self.get_cpx_unsatisfied_cts(self.cts, sol, self.eps)
        for ct, cut, sense, rhs in unsats:
            #Method add() here is CPLEX's CutCallback.add()
            self.add_local(cut, sense, rhs)
            self.nb_cuts += 1
            #print('-- add new cut [{0}]: [{1!s}]'.format(self.nb_cuts, ct))
                
            # to aviod burdensome cuts callback, deleting previous cuts when adding new ones 
            #for i in range(0, int(len(self.cts)*self.frac_pop)):
            #self.cts.pop(0)
          
class LazyCallback(ConstraintCallbackMixin, LazyConstraintCallback):
    #defining lazycallback class for this model
    
    def __init__(self, env):
        #initilize the class
        LazyConstraintCallback.__init__(self, env)
        ConstraintCallbackMixin.__init__(self)
        
        #count the number of lazy cuts added to the model
        self.nb_lazy_cts = 0
        
    def add_lazy_constraints(self, cts):
        self.register_constraints(cts)
        
    @print_called('--> lazy constriant callback called: #{0}')
    def __call__(self):
        #fetch variable values into a solution
        sol = self.make_complete_solution()
        
        #checking the unsatisfied lazy cut and add it to the subproblem
        unsats = self.get_cpx_unsatisfied_cts(self.cts, sol, tolerance = 1e-6)
        for ct, cpx_lhs, sense, cpx_rhs in unsats:
            self.add(cpx_lhs, sense, cpx_rhs)
            self.nb_lazy_cts += 1
            #print(' -- new lazy constraint [{0}]: {1!s}'.format(self.nb_lazy_cts, ct))
'''
class UserCallback(ConstraintCallbackMixin, UserCutCallback):
    '''usercut callback class by doing separation through self-calculation'''
    
    def __init__(self, env):
        UserCutCallback.__init__(self, env)
        ConstraintCallbackMixin.__init__(self)
        #setting the tolerance of violation
        self.eps = 5
        #count the number of cuts added
        self.nb_cuts = 0
        #set limit for the number of cuts added 
        self.cuts_add_limit = 5000
        #count the number of call
        self.nb_call = 0
        #set limit for the number of call
        self.call_limit = 500
        
    #@print_called('--> usercut callback called: #{0}')    
    def __call__(self):
        '''here we first calculating and checking any violation of each constraint based on the 
        incumbent solution. If there exists any we add it to the subproblem.'''
        
        #here we restrict the number of calls and the number of cuts added
        #skip the separation if not at the end of the cut loop
        if self.nb_call > self.call_limit:
            return
        elif self.nb_cuts > self.cuts_add_limit:
            return 
        elif not self.is_after_cut_loop():
            return
        else:
            self.nb_call += 1
            
            #executing separtaion (3a)
            self.separation_single_var_rampup_3a()
            
            #executing separation (3b)
            self.separation_single_var_rampdown_3b()
            
            #executing separation (5a)
            self.separation_single_var_rampup_5a()
            
            #executing separation (5b)
            self.separation_single_var_rampdown_5b()
            
            #executing seperation (6a)
            self.separation_single_var_rampup_6a(1)
            self.separation_single_var_rampup_6a(0)
            
            #executing seperation (6b)
            self.separation_single_var_rampup_6b(1)
            self.separation_single_var_rampup_6b(0)
            
            #executing separation (9a)
            self.separation_double_var_rampup_9a()
            
            #executing separation (9b)
            self.separation_double_var_rampup_9b()
            
            #executing separation (11a)
            self.separation_double_var_rampup_11a()
            
            #executing separation (11b)
            self.separation_double_var_rampup_11b()
            
            #executing separation (12a)
            self.separation_double_var_rampup_12a()
            
            #executing separation (12b)
            self.separation_double_var_rampup_12b()
            
    def separation_single_var_rampup_3a(self):
        '''single variable ramp up eaxct separation for constriant (3a)'''
        #get the current relaxation solution
        sol = self.make_complete_solution()
        #use start and end to denote the accummulation number of each type of generators for current instance
        start_num = 0
        end_num = 0
        #consider 8 types of generators
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #if the number of current type is greater than 0, update the end_num
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                #if the number of generators of current type i for this instance is 0, skip this process
                k_prime = min(uc_min_on[gen_type] - 1, math.floor((uc_gen_ub[gen_type] - \
                              uc_start_ub[gen_type])/uc_rampup_ub[gen_type]), uc_opt_period - 2)
                #use the s_set to store s_i forall i in [1,k], s_i in [0, k_prime] 
                #where y_(t-s_i) - y_(t-s_i-1) > 0
                s_set = np.array([])
                #use t_most_violate to record the most violated y index
                t_most_violate_index = 0
                w_most_violate_value = 0
                #for t = k_prime + 2 we check each point in [t-k_prime+2, t] for the point y_(t-s_i) - 
                #y_(t-s_i-1) > 0 and store the index (t-s_i) in s_set; in the following steps for t in 
                #[t-k_prime+3, T], we only check the last one in s_set in each run. If t-1-k_prime in s_set,
                #we delete it. If y_(t) - y_(t-1) > 0, we add t to s_set
                t = k_prime + 2
                #search for the point y_(t-s_i) - y_(t-s_i-1) > 0 in [t, t-k_prime] 
                for i in range(0, k_prime+1):
                    if sol[self.y_var[start_num+1, t-i]] - sol[self.y_var[start_num+1, t-i-1]] > 0:
                           s_set = np.append(s_set, t-i)
                #calcualte the violation value for t = k_prime + 2
                w = 0
                #if s_set not empty
                if len(s_set) > 0:
                    w = sol[self.x_var[start_num+1, t]] - uc_gen_ub[gen_type]*sol[self.y_var[start_num+1, t]] + \
                        sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - (t-s)*uc_rampup_ub[gen_type])* \
                        (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s-1]]) for s in s_set)
                            
                #if the violation value is greater than 0, we store current t to t_most_violate_index
                if w > w_most_violate_value:
                    t_most_violate_index = t
                    w_most_violate_value = w
                for t in range(k_prime + 3, uc_opt_period + 1):
                    #if s_set is not empty and the last element is t-1-k_prime we delete it
                    if len(s_set) > 0:
                        if s_set[len(s_set)-1] == t - 1 - k_prime:
                            s_set = np.delete(s_set, -1)
                    #if for current t, y_(t) - y_(t-1) > 0 add t to s_set
                    if sol[self.y_var[start_num+1, t]] - sol[self.y_var[start_num+1, t-1]] > 0:
                        s_set = np.insert(s_set, 0, t)
                    #calculate the violation value
                    w = 0
                    if len(s_set) > 0:
                        w = sol[self.x_var[start_num+1, t]] - uc_gen_ub[gen_type]*sol[self.y_var[start_num+1, t]] \
                        + sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - (t-s)*uc_rampup_ub[gen_type])* \
                        (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s-1]]) for s in s_set)
                    if w > w_most_violate_value:
                        t_most_violate_index = t
                        w_most_violate_value = w
                #if the violation value surpass our tolerance, add it as a cut to the subproblem
                if w_most_violate_value > self.eps:
                    s_set = np.array([], dtype = int)
                    #search for the point y_(t-s_i) - y_(t-s_i-1) > 0 in [t, t-k_prime] for the most violated
                    #point index t_most_violate_index
                    for i in range(0, k_prime+1):
                        if sol[self.y_var[start_num+1,t_most_violate_index-i]] - sol[self.y_var[start_num+1, \
                        t_most_violate_index-i-1]] > 0:
                            s_set = np.append(s_set, i)
                    #add the cut to each generator of current type
                    for g in range(start_num+1, end_num+1):
                        '''
                        self.nb_cuts += 1
                        
                        ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t_most_violate_index] - \
                        uc_gen_ub[gen_type]*self.y_var[g, t_most_violate_index] + sum((uc_gen_ub[gen_type] - \
                        uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])*(self.y_var[g, t_most_violate_index-i] - \
                        self.y_var[g, t_most_violate_index-i-1]) for i in s_set) <= 0)
                        self.add(ct_lhs, ct_sense, ct_rhs)
                        '''
                        for t in range(s_set[-1]+2, uc_opt_period+1):
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t] - \
                            uc_gen_ub[gen_type]*self.y_var[g, t] + sum((uc_gen_ub[gen_type] - \
                            uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])*(self.y_var[g, t] - \
                            self.y_var[g, t]) for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                            
                start_num = end_num
                
    def separation_single_var_rampdown_3b(self):
        '''single variable ramp down exact separation for constriant (3b)'''
        #get the current relaxation solution
        sol = self.make_complete_solution()
        #use start and end to denote the accummulation number of each type of generators for current instance
        start_num = 0
        end_num = 0
        #consider 8 types of generators
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #if the number of generators of current type is greater than 0, update the end_num; else move 
                #the next type
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                #calculate k_prime
                k_prime = min(uc_min_off[gen_type] - 1, math.floor((uc_gen_ub[gen_type] - \
                              uc_start_ub[gen_type])/uc_rampup_ub[gen_type]), uc_opt_period - 2)
                #use the s_set to store s_i forall i in [1,k], s_i in [0, k_prime] 
                #where y_(t+s_i) - y_(t+s_i+1) > 0
                s_set = np.array([])
                #use t_most_violate to record the most violated t index and w_most_violate_value for the violation value
                t_most_violate_index = 0
                w_most_violate_value = 0
                #for t = 1 we check each point in [t, t+k_prime] for the point y_(t+s_i) - 
                #y_(t+s_i+1) > 0 and store the index (t+s_i) in s_set; in the following steps for t in 
                #[2, T-s_k-1], we only check the last one in s_set in each run. If t-1 in s_set,
                #we delete it. If y_(t+k_prime) - y_(t+k_prime+1) > 0, we add t+k_prime to s_set
                t = 1
                #search for the point y_(t+s_i) - y_(t+s_i+1) > 0 in [t, t+k_prime] 
                for i in range(0, k_prime+1):
                    if sol[self.y_var[start_num+1, t+i]] - sol[self.y_var[start_num+1, t+i+1]] > 0:
                           s_set = np.append(s_set, t+i)
                #calcualte the violation value for t = 1
                w = 0
                #if s_set not empty
                if len(s_set) > 0:
                    w = sol[self.x_var[start_num+1, t]] - uc_gen_ub[gen_type]*sol[self.y_var[start_num+1, t]] + \
                        sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - (s-t)*uc_rampdown_ub[gen_type])* \
                        (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s+1]]) for s in s_set)
                            
                #if the violation value is greater than 0, we store current t to t_most_violate_index
                if w > w_most_violate_value:
                    t_most_violate_index = t
                    w_most_violate_value = w
                for t in range(2, uc_opt_period - k_prime):
                    #for t from 2 to uc_opt_period - k_prime -1
                    #if s_set is not empty and the first element is t-1 we delete it
                    if len(s_set) > 0:
                        if s_set[0] == t - 1:
                            s_set = np.delete(s_set, 0)
                    #if for current t+k_prime, y_(t+k_prime) - y_(t+k_prime+1) > 0 add t+k_prime to s_set
                    if sol[self.y_var[start_num+1, t+k_prime]] - sol[self.y_var[start_num+1, t+k_prime+1]] > 0:
                        s_set = np.append(s_set, t+k_prime)
                    #calculate the violation value
                    w = 0
                    if len(s_set) > 0:
                        w = sol[self.x_var[start_num+1, t]] - uc_gen_ub[gen_type]*sol[self.y_var[start_num+1, t]] \
                        + sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - (s-t)*uc_rampdown_ub[gen_type])* \
                        (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s+1]]) for s in s_set)
                    if w > w_most_violate_value:
                        t_most_violate_index = t
                        w_most_violate_value = w
                #if the violation value surpass our tolerance, add it as a cut to the subproblem
                if w_most_violate_value > self.eps:
                    s_set = np.array([], dtype = int)
                    #search for the point y_(t+s_i) - y_(t+s_i+1) > 0 in [t, t+k_prime] for the most violated
                    #point index t_most_violate_index
                    for i in range(0, k_prime+1):
                        if sol[self.y_var[start_num+1,t_most_violate_index+i]] - sol[self.y_var[start_num+1, \
                        t_most_violate_index+i+1]] > 0:
                            s_set = np.append(s_set, i)
                    #add the cut to each generator of current type
                    for g in range(start_num+1, end_num+1):
                        '''
                        self.nb_cuts += 1
                        ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t_most_violate_index] - \
                        uc_gen_ub[gen_type]*self.y_var[g,t_most_violate_index] + sum((uc_gen_ub[gen_type] - \
                        uc_start_ub[gen_type] - i*uc_rampdown_ub[gen_type])*(self.y_var[g, t_most_violate_index+i] - \
                        self.y_var[g, t_most_violate_index+i+1]) for i in s_set) <= 0)
                        self.add(ct_lhs, ct_sense, ct_rhs)
                        '''
                        for t in range(1, uc_opt_period - s_set[-1]):
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t] - \
                            uc_gen_ub[gen_type]*self.y_var[g,t] + sum((uc_gen_ub[gen_type] - \
                            uc_start_ub[gen_type] - i*uc_rampdown_ub[gen_type])*(self.y_var[g, t+i] - \
                            self.y_var[g, t+i+1]) for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                start_num = end_num
               
    def separation_single_var_rampup_5a(self):
        '''single variable ramp up exact separation for constraint (5a)'''
        #getting the current relaxation solution
        sol = self.make_complete_solution()
        #use the start and end index to denote the start and end index of current type of generator in current instances
        start_num = 0
        end_num = 0
        #consider the 8 types of generators
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                
                #if the number of current type of generator is not 0, update end_index
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                #calculate beta
                
                beta = min(uc_min_on[gen_type]-1, (uc_gen_ub[gen_type] - \
                                                   uc_start_ub[gen_type])/uc_rampup_ub[gen_type])
                gamma = min(uc_min_on[gen_type]-1, math.floor((uc_gen_ub[gen_type] - \
                       uc_start_ub[gen_type])/uc_rampup_ub[gen_type]), uc_opt_period - 3)
                #set the s_set to store all s_i, i in [1,k], s_i in [0,beta], satisfying y_var[t-s_i] - 
                #y_var[t-s_i-1] > 0
                s_set = np.array([])
                
                #use t_most_violate_index and w_most_violate_value to store the most violated index of t and it 
                #violation value
                t_most_violate_index = 0
                w_most_violate_value = 0
                #for t = beta + 2, check for s_i, y_var[t-s_i] - y_var[t-s_i-1] > 0, and store t-s_i in s_set
                t = gamma + 2
                #search for t-s_i such that y_var[t-s_i] - y_var[t-s_i-1] > 0, for s_i on [0, beta]
                for i in range(0, gamma+1):
                    if sol[self.y_var[start_num+1, t-i]] - sol[self.y_var[start_num+1, t-i-1]] > 0:
                        s_set = np.append(s_set, t-i)
                #calculate the violation value for t = gamma + 2
                
                w = 0
                if len(s_set) > 0:
                    #if s_set is not empty
                    w = sol[self.x_var[start_num+1, t]] - (uc_gen_ub[gen_type] - beta*uc_rampup_ub[gen_type])* \
                        sol[self.y_var[start_num+1, t]] - beta*uc_rampup_ub[gen_type]*sol[self.y_var[start_num+1, \
                        t+1]] + sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - (t-s)*uc_rampup_ub[gen_type])* \
                        (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s-1]]) for s in s_set)
                if w > w_most_violate_value:
                    t_most_violate_index = t
                    w_most_violate_value = w
                
                #For t = gamma+3,...,T-1, we only check the first element in s_set
                #if it is t-1-gamma, delete it from s_set; And then check whether y_var[t] - y_var[t-1] > 0, if yes
                #add t into s_set
                for t in range(gamma+3, uc_opt_period):
                    #if s_set in not empty, check whether the last element is t-1-beta, if so, delete it
                    if len(s_set) > 0:
                        if s_set[len(s_set)-1] == t - 1 - gamma:
                            s_set = np.delete(s_set, -1)
                    if sol[self.y_var[start_num+1, t]] - sol[self.y_var[start_num+1, t-1]] > 0:
                    #if for the currrent t, y_var[t] - y_var[t-1] > 0, insert t at the first position of s_set
                        s_set = np.insert(s_set, 0, t)
                    #calculate the violation value
                    w = 0
                    if len(s_set) > 0:
                        w = sol[self.x_var[start_num+1, t]] - (uc_gen_ub[gen_type] - beta*uc_rampup_ub[gen_type])* \
                            sol[self.y_var[start_num+1, t]] - beta*uc_rampup_ub[gen_type]*sol[self.y_var[start_num+1, \
                            t+1]] + sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - (t-s)*uc_rampup_ub[gen_type])* \
                            (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s-1]]) for s in s_set)
                    if w > w_most_violate_value:
                        #if the current violation is better than w_most_violate_value, update most violate
                        t_most_violate_index = t
                        w_most_violate_value = w
                
                if w_most_violate_value > self.eps:
                    #if the violation exceeds our tolerance range, find s_set for t_most_violate_index, 
                    #forall s_i in [0,gamma] satisfying, y_var[t_most_violate_index - s_i] - 
                    #y_var[t_most_violate_index - s_i -1] > 0
                    s_set = np.array([], dtype = int)
                    for i in range(0, gamma+1):
                        if sol[self.y_var[start_num+1, t_most_violate_index-i]] - sol[self.y_var[start_num+1, \
                        t_most_violate_index-i-1]] > 0:
                            s_set = np.append(s_set, i)
                    #add this cut to all generators of current type
                    for g in range(start_num+1, end_num+1):
                        '''
                        self.nb_cuts += 1
                        ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t_most_violate_index] - \
                        (uc_gen_ub[gen_type] - beta*uc_rampup_ub[gen_type])*self.y_var[g, t_most_violate_index] - \
                        beta*uc_rampup_ub[gen_type]*self.y_var[g, t_most_violate_index+1] + \
                        sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])* \
                        (self.y_var[g, t_most_violate_index-i] - self.y_var[g, t_most_violate_index-i-1]) \
                        for i in s_set) <= 0)
                        self.add(ct_lhs, ct_sense, ct_rhs)
                        '''
                        for t in range(s_set[-1]+2, uc_opt_period):
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t] - \
                            (uc_gen_ub[gen_type] - beta*uc_rampup_ub[gen_type])*self.y_var[g, t] - \
                            beta*uc_rampup_ub[gen_type]*self.y_var[g, t+1] + \
                            sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])* \
                            (self.y_var[g, t-i] - self.y_var[g, t-i-1]) \
                            for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                start_num = end_num
                
    def separation_single_var_rampdown_5b(self):
        '''single variable ramp down exact separation for constraint (5b)'''
        #getting the current relaxation solution
        sol = self.make_complete_solution()
        #use the start and end index to denote the start and end index of current type of generator in current instances
        start_num = 0
        end_num = 0
        #consider the 8 types of generators
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #if the number of current type of generator is not 0, update end_index
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                #calculate beta and gamma
                beta = min(uc_min_on[gen_type]-1, (uc_gen_ub[gen_type] - \
                                                   uc_start_ub[gen_type])/uc_rampup_ub[gen_type])
                gamma = min(uc_min_on[gen_type]-1, math.floor((uc_gen_ub[gen_type] - \
                       uc_start_ub[gen_type])/uc_rampup_ub[gen_type]), uc_opt_period - 3)
                #set the s_set to store all s_i, i in [1,k], s_i in [0,gamma], satisfying y_var[t+s_i] - 
                #y_var[t+s_i+1] > 0
                s_set = np.array([])
                #use t_most_violate_index and w_most_violate_value to store the most violated index of t and it 
                #violation value
                t_most_violate_index = 0
                w_most_violate_value = 0
                #for t = 2, check for s_i, y_var[t+s_i] - y_var[t+s_i+1] > 0, and store t+s_i in s_set
                t = 2
                #search for t-s_i such that y_var[t+s_i] - y_var[t+s_i+1] > 0, for s_i on [0, beta]
                for i in range(0, gamma+1):
                    if sol[self.y_var[start_num+1, t+i]] - sol[self.y_var[start_num+1, t+i+1]] > 0:
                        s_set = np.append(s_set, t+i)
                #calculate the violation value for t = 2
                w = 0
                if len(s_set) > 0:
                    #if s_set is not empty
                    w = sol[self.x_var[start_num+1, t]] - (uc_gen_ub[gen_type] - beta*uc_rampup_ub[gen_type])* \
                        sol[self.y_var[start_num+1, t]] - beta*uc_rampup_ub[gen_type]*sol[self.y_var[start_num+1, \
                        t-1]] + sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - (s-t)*uc_rampup_ub[gen_type])* \
                        (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s+1]]) for s in s_set)
                if w > w_most_violate_value:
                    t_most_violate_index = t
                    w_most_violate_value = w
                #For t = 3,...,T-gamma-1, we only check the first element in s_set
                #if it is t-1-gamma, delete it from s_set; And then check whether y_var[t] - y_var[t-1] > 0, if yes
                #add t into s_set
                for t in range(3, uc_opt_period-gamma):
                    #if s_set in not empty, check whether the last element is t-1-gamma, if so, delete it
                    if len(s_set) > 0:
                        if s_set[0] == t - 1:
                            s_set = np.delete(s_set, 0)
                    if sol[self.y_var[start_num+1, t+gamma]] - sol[self.y_var[start_num+1, t+gamma+1]] > 0:
                    #if for the currrent t, y_var[t] - y_var[t-1] > 0, insert t at the first position of s_set
                        s_set = np.append(s_set, t+gamma)
                    #calculate the violation value
                    w = 0
                    if len(s_set) > 0:
                        w = sol[self.x_var[start_num+1, t]] - (uc_gen_ub[gen_type] - beta*uc_rampup_ub[gen_type])* \
                            sol[self.y_var[start_num+1, t]] - beta*uc_rampup_ub[gen_type]*sol[self.y_var[start_num+1, \
                            t-1]] + sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - (s-t)*uc_rampup_ub[gen_type])* \
                            (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s+1]]) for s in s_set)
                    if w > w_most_violate_value:
                        #if the current violation is better than w_most_violate_value, update most violate
                        t_most_violate_index = t
                        w_most_violate_value = w
                if w_most_violate_value > self.eps:
                    #if the violation exceeds our tolerance range, find s_set for t_most_violate_index, 
                    #forall s_i in [0,beta] satisfying, y_var[t_most_violate_index - s_i] - 
                    #y_var[t_most_violate_index - s_i -1] > 0
                    s_set = np.array([], dtype = int)
                    for i in range(0, gamma+1):
                        if sol[self.y_var[start_num+1, t_most_violate_index+i]] - sol[self.y_var[start_num+1, \
                        t_most_violate_index+i+1]] > 0:
                            s_set = np.append(s_set, i)
                    #add this cut to all generators of current type
                    for g in range(start_num+1, end_num+1):
                        '''
                        self.nb_cuts += 1
                        ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t_most_violate_index] - \
                        (uc_gen_ub[gen_type] - beta*uc_rampup_ub[gen_type])*self.y_var[g, t_most_violate_index] -\
                        beta*uc_rampup_ub[gen_type]*self.y_var[g, t_most_violate_index-1] + \
                        sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])* \
                        (self.y_var[g, t_most_violate_index+i] - self.y_var[g, t_most_violate_index+i+1]) \
                        for i in s_set) <= 0)
                        self.add(ct_lhs, ct_sense, ct_rhs)
                        '''
                        for t in range(2, uc_opt_period - int(s_set[-1])):
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t] - \
                            (uc_gen_ub[gen_type] - beta*uc_rampup_ub[gen_type])*self.y_var[g, t] -\
                            beta*uc_rampup_ub[gen_type]*self.y_var[g, t-1] + \
                            sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - i*uc_rampup_ub[gen_type])* \
                            (self.y_var[g, t+i] - self.y_var[g, t+i+1]) \
                            for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                start_num = end_num
    
    def separation_single_var_rampup_6a(self, token):
        '''single variable ramp up exact separation for constraint (6a)'''
        #getting the current relaxation solution
        sol = self.make_complete_solution()
        #use the start and end index to denote the start and end index of current type of generator in current instances
        start_num = 0
        end_num = 0
        #consider the 8 types of generators
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #if the number for generators of current type is nonzero, update end_index
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                #set the value of beta based on token: if token is 1, set beta = min(L, (Cupper - Vupper)/V)
                #if token is 0, set beta = 0. Both are facet defining for conv(P)
                if token == 1:
                    beta = min(uc_min_on[gen_type], (uc_gen_ub[gen_type] - uc_start_ub[gen_type])/uc_rampup_ub[gen_type])
                elif token == 0:
                    beta = 0
                else:
                    raise ValueError('Wrong token was given for callback constraints 6a')
                #set the value of gamma
                gamma = min(uc_min_on[gen_type], math.floor((uc_gen_ub[gen_type] - \
                        uc_start_ub[gen_type])/uc_rampup_ub[gen_type]), uc_opt_period-2)
                #set the set S to store all t-i, i in [1, gamma] such that y_[t-i] - y_[t-i-1] > 0 for current t
                s_set = np.array([])
                #we use t_most_violate_index and value to store the index and violation value of the 
                #greatest violation of (10a)
                t_most_violate_index = 0
                w_most_violate_value = 0
                #for t = gamma+2, we find all t-i, i in [1,gamma], such that y_[t-i] - y_[t-i-1] > 0 and store
                #them in s_set
                t = gamma + 2
                for i in range(1, gamma+1):
                    if sol[self.y_var[start_num+1,t-i]] - sol[self.y_var[start_num+1, t-i-1]] > 0:
                        s_set = np.append(s_set, t-i)
                #calculate the violation value for t = gamma+2
                w = 0
                if len(s_set) > 0:
                    # if the s_set is not empty
                    w = sol[self.x_var[start_num+1, t]] - (uc_start_ub[gen_type]+beta*uc_rampup_ub[gen_type])*\
                        sol[self.y_var[start_num+1, t]] - (uc_gen_ub[gen_type]-uc_start_ub[gen_type] - \
                        beta*uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t-1]] + sum((uc_gen_ub[gen_type] \
                        - uc_start_ub[gen_type] - (t-s)*uc_rampup_ub[gen_type])*(sol[self.y_var[start_num+1,s]] - \
                        sol[self.y_var[start_num+1, s-1]]) for s in s_set)
                if w > w_most_violate_value:
                    t_most_violate_index = t
                    w_most_violate_value = w
                #for t = gamma+3,...,T, we only check whether y_[t-1] - y_[t-2] > 0. If it is we add t-1 to the
                #first place of s_set; and then check whether t-1-gamma is in s_set, if so we delete it from s_set
                for t in range(gamma+3, uc_opt_period+1):
                    #for all t =gamma+3,...,T
                    if len(s_set) > 0:
                        #if s_set is not emoty, check whether its last element is t-1-gamma, if yes, delete it
                        if s_set[len(s_set)-1] == t-1-gamma:
                            s_set = np.delete(s_set, -1)
                    if sol[self.y_var[start_num+1, t-1]] - sol[self.y_var[start_num+1, t-2]] > 0:
                        #if y_[t-1] - y_[t-2] > 0, insert t-1 to the first place of s_set
                        s_set = np.insert(s_set, 0, t-1)
                    #calculate the violation value for current t
                    w = 0
                    if len(s_set) > 0:
                        # if the s_set is not empty
                        w = sol[self.x_var[start_num+1, t]] - (uc_start_ub[gen_type]+beta*uc_rampup_ub[gen_type])*\
                            sol[self.y_var[start_num+1, t]] - (uc_gen_ub[gen_type]-uc_start_ub[gen_type] - \
                            beta*uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t-1]] + sum((uc_gen_ub[gen_type] \
                            - uc_start_ub[gen_type] - (t-s)*uc_rampup_ub[gen_type])*(sol[self.y_var[start_num+1,s]] - \
                            sol[self.y_var[start_num+1, s-1]]) for s in s_set)
                    if w > w_most_violate_value:
                        t_most_violate_index = t
                        w_most_violate_value = w
                if w_most_violate_value > self.eps:
                    #if the most violation exceeds our tolerance add it to all generators of current type
                    s_set = np.array([], dtype = int)
                    for i in range(1, gamma+1):
                        if sol[self.y_var[start_num+1, t_most_violate_index-i]] - sol[self.y_var[start_num+1, \
                                t_most_violate_index-i-1]] > 0:
                            s_set = np.append(s_set, i)
                    #add this constraint to all generators of current type
                    for g in range(start_num+1, end_num+1):
                        '''
                        self.nb_cuts += 1
                        #updating the counting of number of cuts we added
                        ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g,t_most_violate_index] - \
                            (uc_start_ub[gen_type] + beta*uc_rampup_ub[gen_type])*self.y_var[g, t_most_violate_index] - \
                            (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - beta*uc_rampup_ub[gen_type])* \
                            self.y_var[g,t_most_violate_index-1] + sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - \
                            i*uc_rampup_ub[gen_type])*(self.y_var[g, t_most_violate_index-i] - \
                            self.y_var[g, t_most_violate_index-i-1]) for i in s_set) <= 0)
                        self.add(ct_lhs, ct_sense, ct_rhs)
                        '''
                        for t in range(int(s_set[-1])+2, uc_opt_period+1):
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g,t] - \
                                (uc_start_ub[gen_type] + beta*uc_rampup_ub[gen_type])*self.y_var[g, t] - \
                                (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - beta*uc_rampup_ub[gen_type])* \
                                self.y_var[g,t-1] + sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - \
                                i*uc_rampup_ub[gen_type])*(self.y_var[g, t-i] - \
                                self.y_var[g, t-i-1]) for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                #move the start index to the end index of current type of generators
                start_num = end_num
                
    def separation_single_var_rampup_6b(self, token):
        '''single variable ramp up exact separation for constraint (6b)'''
        #getting the current relaxation solution
        sol = self.make_complete_solution()
        #use the start and end index to denote the start and end index of current type of generator in current instances
        start_num = 0
        end_num = 0
        #consider the 8 types of generators
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #if the number for generators of current type is nonzero, update end_index
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                #set the value of beta based on token: if token is 1, set beta = min(L, (Cupper - Vupper)/V)
                #if token is 0, set beta = 0. Both are facet defining for conv(P)
                if token == 1:
                    beta = min(uc_min_on[gen_type], (uc_gen_ub[gen_type] - \
                                                     uc_start_ub[gen_type])/uc_rampup_ub[gen_type])
                elif token == 0:
                    beta = 0
                else:
                    return
                #set the value of gamma
                gamma = min(uc_min_on[gen_type], math.floor((uc_gen_ub[gen_type] - \
                        uc_start_ub[gen_type])/uc_rampup_ub[gen_type]), uc_opt_period-2)
                #set the set S to store all t-i, i in [1, gamma] such that y_[t+i] - y_[t+i+1] > 0 for current t
                s_set = np.array([])
                #we use t_most_violate_index and value to store the index and violation value of the 
                #greatest violation of (10a)
                t_most_violate_index = 0
                w_most_violate_value = 0
                #for t = gamma+2, we find all t-i, i in [1,gamma], such that y_[t+i] - y_[t+i+1] > 0 and store
                #them in s_set
                t = 1
                for i in range(1, gamma+1):
                    if sol[self.y_var[start_num+1,t+i]] - sol[self.y_var[start_num+1, t+i+1]] > 0:
                        s_set = np.append(s_set, t+i)
                #calculate the violation value for t = 1
                w = 0
                if len(s_set) > 0:
                    # if the s_set is not empty
                    w = sol[self.x_var[start_num+1, t]] - (uc_start_ub[gen_type]+beta*uc_rampup_ub[gen_type])*\
                        sol[self.y_var[start_num+1, t]] - (uc_gen_ub[gen_type]-uc_start_ub[gen_type] - \
                        beta*uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t+1]] + sum((uc_gen_ub[gen_type] \
                        - uc_start_ub[gen_type] - (s-t)*uc_rampup_ub[gen_type])*(sol[self.y_var[start_num+1,s]] - \
                        sol[self.y_var[start_num+1, s+1]]) for s in s_set)
                if w > w_most_violate_value:
                    t_most_violate_index = t
                    w_most_violate_value = w
                #for each t = 2,...,T-gamma-1, we only check whether y_[t+1+gamma] - y_[t+1+gamma+1] > 0. If it is we add t+1+gamma 
                #to the last place of s_set; and then check whether t+1 is in s_set, if so we delete it from s_set
                for t in range(2, uc_opt_period-gamma):
                    #for all t =2,...,T-gamma-1
                    if len(s_set) > 0:
                        #if s_set is not emoty, check whether its first element is t, if yes, delete it
                        if s_set[0] == t:
                            s_set = np.delete(s_set, 0)
                    if sol[self.y_var[start_num+1, t+gamma]] - sol[self.y_var[start_num+1, t+gamma+1]] > 0:
                        #if y_[t+gamma] - y_[t+gamma] > 0, insert t-1 to the first place of s_set
                        s_set = np.append(s_set, t+gamma)
                    #calculate the violation value for current t
                    w = 0
                    if len(s_set) > 0:
                        # if the s_set is not empty
                        w = sol[self.x_var[start_num+1, t]] - (uc_start_ub[gen_type]+beta*uc_rampup_ub[gen_type])*\
                            sol[self.y_var[start_num+1, t]] - (uc_gen_ub[gen_type]-uc_start_ub[gen_type]- \
                            beta*uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t+1]] + sum((uc_gen_ub[gen_type] \
                            - uc_start_ub[gen_type] - (s-t)*uc_rampup_ub[gen_type])*(sol[self.y_var[start_num+1,s]] - \
                            sol[self.y_var[start_num+1, s+1]]) for s in s_set)
                    if w > w_most_violate_value:
                        t_most_violate_index = t
                        w_most_violate_value = w
                if w_most_violate_value > self.eps:
                    #if the most violation exceeds our tolerance add it to all generators of current type
                    s_set = np.array([], dtype = int)
                    for i in range(1, gamma+1):
                        if sol[self.y_var[start_num+1, t_most_violate_index+i]] - sol[self.y_var[start_num+1, \
                                t_most_violate_index+i+1]] > 0:
                            s_set = np.append(s_set, i)
                    #add this constraint to all generators of current type
                    for g in range(start_num+1, end_num+1):
                        '''
                        self.nb_cuts += 1
                        #updating the counting of number of cuts we added
                        ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g,t_most_violate_index] - \
                            (uc_start_ub[gen_type] + beta*uc_rampup_ub[gen_type])*self.y_var[g, t_most_violate_index] - \
                            (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - beta*uc_rampup_ub[gen_type])* \
                            self.y_var[g,t_most_violate_index+1] + sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - \
                            i*uc_rampup_ub[gen_type])*(self.y_var[g, t_most_violate_index+i] - \
                            self.y_var[g, t_most_violate_index+i+1]) for i in s_set) <= 0)
                        self.add(ct_lhs, ct_sense, ct_rhs)
                        '''
                        for t in range(1, uc_opt_period-s_set[-1]):
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g,t] - \
                                (uc_start_ub[gen_type] + beta*uc_rampup_ub[gen_type])*self.y_var[g, t] - \
                                (uc_gen_ub[gen_type] - uc_start_ub[gen_type] - beta*uc_rampup_ub[gen_type])* \
                                self.y_var[g,t+1] + sum((uc_gen_ub[gen_type] - uc_start_ub[gen_type] - \
                                i*uc_rampup_ub[gen_type])*(self.y_var[g, t+i] - \
                                self.y_var[g, t+i+1]) for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                #move the start index to the end index of current type of generators
                start_num = end_num
                        
    def separation_double_var_rampup_9a(self):
        '''double var ramp up upper bound exact separation for constraint (9a)'''
        #getting the current relaxation solution
        sol = self.make_complete_solution()
        #use the start_num and end_num to denote the start and end index for the current type of generator in current instance
        start_num = 0
        end_num = 0
        #consider the 8 types of generators
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #if the number of current type of generators in current instance is not 0, update end index
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                #set the values of all parameters
                #use w and t to record the largest violation of 9a among all T periods for current type of generators
                w_most_violate_value = 0
                t_most_violate_index = 0
                k_most_violate_index = 0
                for t in range(2, uc_opt_period+1):
                    # set the value of k'
                    k_prime = min(t-1, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampup_ub[gen_type]))
                    #use w_hat, k_hat, and s_set to store the largest violation value, index, and S set for current 
                    #type of generator under current t
                    w_hat = 0
                    k_hat = 0
                    s_set = np.array([])
                    w = 0
                    #use temp to denote the temporary value of x_t - x_[t-k] - (uc_gen_lb + kV)*y_t + uc_gen_lb*y_[t-k] to find k0
                    #k0 is the start value of k
                    k0 = 1
                    temp0 =  sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t-k0]] - (uc_gen_lb[gen_type] + \
                            k0*uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t]] + \
                            uc_gen_lb[gen_type]*sol[self.y_var[start_num+1, t-k0]]
                    for i in range(1, k_prime+1):
                        temp1 = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t-i]] - (uc_gen_lb[gen_type] + \
                                i*uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t]] + \
                                uc_gen_lb[gen_type]*sol[self.y_var[start_num+1, t-i]]
                        if temp1 > temp0:
                            #update k0
                            k0 = i
                            temp0 = temp1
                    for i in range(0, min(k0, uc_min_on[gen_type])):
                        #find i between 0 and min(k0-1, L-1) such that y_[t-i] - y_[t-i-1] > 0
                        if sol[self.y_var[start_num+1, t-i]] - sol[self.y_var[start_num+1, t-i-1]] > 0:
                            s_set = np.append(s_set, t-i)
                    if len(s_set) > 0:
                        w = temp0 + sum((uc_gen_lb[gen_type] + (k0 - t + s)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*\
                            (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s-1]]) for s in s_set)
                    if w > w_hat:
                        # update w_hat
                        w_hat = w
                        k_hat = k0
                    if k0 < k_prime:
                        #if k0 is not equal to the largest possible value of k, check for violation from k0+1 to k_prime
                        for i in range(k0+1, k_prime+1):
                            w = 0
                            if i <= uc_min_on[gen_type] and sol[self.y_var[start_num+1, t-i+1]] - sol[self.y_var[start_num+1, t-i]] > 0:
                                #if i-1 <= min(L-1, k-1), add it to s_set. Here k = i
                                s_set = np.append(s_set, t-i+1)
                            if len(s_set) > 0:
                                w = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t-i]] - (uc_gen_lb[gen_type] + \
                                    i*uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t]] + uc_gen_lb[gen_type]*sol[self.y_var[start_num+1, \
                                    t-i]] + sum((uc_gen_lb[gen_type] + (i-t+s)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*\
                                    (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s-1]]) for s in s_set)
                            if w > w_hat:
                                #update w_hat
                                w_hat = w
                                k_hat = i
                    if w_hat > w_most_violate_value:
                        #update the largest violation value and index
                        w_most_violate_value = w_hat
                        t_most_violate_index = t
                        k_most_violate_index = k_hat
                if w_most_violate_value > self.eps:
                    s_set = np.array([], dtype = int)
                    for i in range(0, min(k_most_violate_index, uc_min_on[gen_type])):
                        if sol[self.y_var[start_num+1, t_most_violate_index-i]] - sol[self.y_var[start_num+1, t_most_violate_index-i-1]] > 0:
                            s_set = np.append(s_set, i)
                    for g in range(start_num+1, end_num+1):
                        '''
                        #add the most violation constraint to all current type of generators
                        self.nb_cuts += 1
                        #update the number of cuts we added
                        ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t_most_violate_index] -\
                            self.x_var[g, t_most_violate_index-k_most_violate_index] - (uc_gen_lb[gen_type] + k_most_violate_index*\
                            uc_rampup_ub[gen_type])*self.y_var[g, t_most_violate_index] + uc_gen_lb[gen_type]*self.y_var[g, \
                            t_most_violate_index - k_most_violate_index] + sum((uc_gen_lb[gen_type]+(k_most_violate_index-i)*\
                            uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*(self.y_var[g, t_most_violate_index-i] - self.y_var[g, \
                            t_most_violate_index-i-1]) for i in s_set) <= 0)
                        self.add(ct_lhs, ct_sense, ct_rhs)
                        '''
                        for t in range(k_most_violate_index+1, uc_opt_period+1):
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t] -\
                                self.x_var[g, t-k_most_violate_index] - (uc_gen_lb[gen_type] + k_most_violate_index*\
                                uc_rampup_ub[gen_type])*self.y_var[g, t] + uc_gen_lb[gen_type]*self.y_var[g, \
                                t - k_most_violate_index] + sum((uc_gen_lb[gen_type]+(k_most_violate_index-i)*\
                                uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*(self.y_var[g, t-i] - self.y_var[g, \
                                t-i-1]) for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                #move the start index to the end index of current type of generators
                start_num = end_num
        
    def separation_double_var_rampup_9b(self):
        '''double var ramp down exact separation for constraint (9b)'''
        #getting the current relaxation solution
        sol = self.make_complete_solution()
        #use the start_num and end_num to denote the start and end index for the current type of generator in current instance
        start_num = 0
        end_num = 0
        #consider the 8 types of generators
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #if the number of current type of generators in current instance is not 0, update end index
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                #set the values of all parameters
                #use w and t to record the largest violation of 9a among all T periods for current type of generators
                w_most_violate_value = 0
                t_most_violate_index = 0
                k_most_violate_index = 0
                for t in range(1, uc_opt_period-1):
                    #set the value of k'
                    k_prime = min(uc_opt_period-t, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampdown_ub[gen_type]))
                    #use w_hat, k_hat, and s_set to store the largest violation value, index, and S set for current 
                    #type of generator under current t
                    w_hat = 0
                    k_hat = 0
                    s_set = np.array([])
                    w = 0
                    #use temp to denote the temporary value of x_t - x_[t+k] - (uc_gen_lb + kV)*y_t + uc_gen_lb*y_[t+k] to find k0
                    #k0 is the start value of k
                    k0 = 1
                    temp0 =  sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t+k0]] - (uc_gen_lb[gen_type] + \
                            k0*uc_rampdown_ub[gen_type])*sol[self.y_var[start_num+1, t]] + \
                            uc_gen_lb[gen_type]*sol[self.y_var[start_num+1, t+k0]]
                    for i in range(1, k_prime+1):
                        temp1 = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t+i]] - (uc_gen_lb[gen_type] + \
                                i*uc_rampdown_ub[gen_type])*sol[self.y_var[start_num+1, t]] + \
                                uc_gen_lb[gen_type]*sol[self.y_var[start_num+1, t+i]]
                        if temp1 > temp0:
                            #update k0
                            k0 = i
                            temp0 = temp1
                    for i in range(0, min(k0, uc_min_on[gen_type])):
                        #find i between 0 and min(k0-1, L-1) such that y_[t+i] - y_[t+i+1] > 0
                        if sol[self.y_var[start_num+1, t+i]] - sol[self.y_var[start_num+1, t+i+1]] > 0:
                            s_set = np.append(s_set, t+i)
                    if len(s_set) > 0:
                        w = temp0 + sum((uc_gen_lb[gen_type] + (k0 - s + t)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*\
                            (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s+1]]) for s in s_set)
                    if w > w_hat:
                        # update w_hat
                        w_hat = w
                        k_hat = k0
                    if k0 < k_prime:
                        #if k0 is not equal to the largest possible value of k, check for violation from k0+1 to k_prime
                        for i in range(k0+1, k_prime+1):
                            w = 0
                            if i <= uc_min_on[gen_type] and sol[self.y_var[start_num+1, t+i-1]] - sol[self.y_var[start_num+1, t+i]] > 0:
                                #if i-1 <= min(L-1, k-1), add it to s_set. Here k = i
                                s_set = np.append(s_set, t+i-1)
                            if len(s_set) > 0:
                                w = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t+i]] - (uc_gen_lb[gen_type] + \
                                    i*uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t]] + uc_gen_lb[gen_type]*sol[self.y_var[start_num+1, \
                                    t+i]] + sum((uc_gen_lb[gen_type] + (i-s+t)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*\
                                    (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s+1]]) for s in s_set)
                            if w > w_hat:
                                #update w_hat
                                w_hat = w
                                k_hat = i
                    if w_hat > w_most_violate_value:
                        #update the largest violation value and index
                        w_most_violate_value = w_hat
                        t_most_violate_index = t
                        k_most_violate_index = k_hat
                if w_most_violate_value > self.eps:
                    s_set = np.array([], dtype = int)
                    for i in range(0, min(k_most_violate_index, uc_min_on[gen_type])):
                        if sol[self.y_var[start_num+1, t_most_violate_index+i]] - sol[self.y_var[start_num+1, t_most_violate_index+i+1]] > 0:
                            s_set = np.append(s_set, i)
                    for g in range(start_num+1, end_num+1):
                        '''
                        #add the most violation constraint to all current type of generators
                        self.nb_cuts += 1
                        #update the number of cuts we added
                        ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t_most_violate_index] -\
                            self.x_var[g, t_most_violate_index + k_most_violate_index] - (uc_gen_lb[gen_type] + k_most_violate_index*\
                            uc_rampdown_ub[gen_type])*self.y_var[g, t_most_violate_index] + uc_gen_lb[gen_type]*self.y_var[g, \
                            t_most_violate_index + k_most_violate_index] + sum((uc_gen_lb[gen_type]+(k_most_violate_index-i)*\
                            uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])*(self.y_var[g, t_most_violate_index+i] - self.y_var[g, \
                            t_most_violate_index+i+1]) for i in s_set) <= 0)
                        self.add(ct_lhs, ct_sense, ct_rhs)
                        '''
                        for t in range(1, uc_opt_period-k_most_violate_index+1):
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t] -\
                                self.x_var[g, t + k_most_violate_index] - (uc_gen_lb[gen_type] + k_most_violate_index*\
                                uc_rampdown_ub[gen_type])*self.y_var[g, t] + uc_gen_lb[gen_type]*self.y_var[g, \
                                t + k_most_violate_index] + sum((uc_gen_lb[gen_type]+(k_most_violate_index-i)*\
                                uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])*(self.y_var[g, t+i] - self.y_var[g, \
                                t+i+1]) for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                #move the start index to the end index of current type of generators
                start_num = end_num 
                        
    def separation_double_var_rampup_11a(self):
        '''double var exact separation for constraint (11a)'''
        #get the current relaxation solution
        sol = self.make_complete_solution()
        #use start_num and end_num to denote the start index and end index of current type of generators under 
        #current instance
        start_num = 0
        end_num = 0
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #the number of current type of generators for current instance should be greater than 0
                #update the end_num
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                if uc_min_on[gen_type] > 1:
                    # the minimum on period should be greater than 1, i.e., L > 1
                    #define the most violation value, indexes of current relaxation solution
                    w_most_violate_value = 0
                    t_most_violate_index = 0
                    k_most_violate_index = 0
                    m_most_violate_index = 0
                    #calculate the largest violation for each t in [2, T-1]
                    for t in range(2, uc_opt_period):
                        #set the value of k', i.e., the largest possible value of k
                        k_prime = min(t-1, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampup_ub[gen_type]))
                        #define the largest violation value and index for given t
                        w_hat = 0
                        k_hat = 0
                        m_hat = 0
                        for k in range(1, k_prime+1):
                            #define the upper bound of s_q, i.e., alpha, and the upper bound of m, i.e., beta
                            alpha = min(k-1, uc_min_on[gen_type]-2)
                            beta = min(alpha, uc_opt_period - t -1)
                            #define the largest violation value and indexes for given t and given k
                            w_bar = 0
                            m_bar = 0
                            s_set = np.array([])
                            #find the value of m0
                            m0 = 0
                            temp0 = -(uc_gen_lb[gen_type] + (k-m0)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*sol[self.y_var[start_num+1, \
                                    t+m0+1]] - uc_rampup_ub[gen_type]*sum(sol[self.y_var[start_num+1, t+j]] for j in range(1, m0+1))
                            for m in range(0, beta+1):
                                temp1 = -(uc_gen_lb[gen_type] + (k-m)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])* \
                                      sol[self.y_var[start_num+1, t+m+1]] - uc_rampup_ub[gen_type]*sum(sol[self.y_var[start_num+1, t+j]] for \
                                      j in range(1, m+1))
                                if temp1 > temp0:
                                    temp0 = temp1
                                    m0 = m
                            for i in range(0, alpha-m0+1):
                                #search for i between 0 and alpha - m0 such that y_[t-i] - y_[t-i-1] > 0
                                if sol[self.y_var[start_num+1, t-i]] - sol[self.y_var[start_num+1, t-i-1]] > 0:
                                    s_set = np.append(s_set, t-i)
                            #calculate the violation value for m0
                            w = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t-k]] - (uc_gen_lb[gen_type] + \
                                (k-m0)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*sol[self.y_var[start_num+1, t+m0+1]] - \
                                uc_rampup_ub[gen_type]*sum(sol[self.y_var[start_num+1, t+j]] for j in range(1, m0+1)) - \
                                uc_start_ub[gen_type]*sol[self.y_var[start_num+1, t]] + uc_gen_lb[gen_type]*sol[self.y_var[start_num+1, \
                                t-k]] + sum((uc_gen_lb[gen_type] + (k-t+s)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])* \
                                (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s-1]]) for s in s_set)
                            if w > w_bar:
                                w_bar = w
                                m_bar = m0
                            if m0 > 0:
                                #if m0 is not 0, we check for any better violation of m in [0,m-1]
                                m = m0-1
                                while m >= 0:
                                    if sol[self.y_var[start_num+1, t-alpha+m]] - sol[self.y_var[start_num+1, t-alpha+m-1]] > 0:
                                        s_set = np.append(s_set, t-alpha+m)
                                    w = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t-k]] - (uc_gen_lb[gen_type] + \
                                        (k-m)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*sol[self.y_var[start_num+1, t+m+1]] - \
                                        uc_rampup_ub[gen_type]*sum(sol[self.y_var[start_num+1, t+j]] for j in range(1, m+1)) - \
                                        uc_start_ub[gen_type]*sol[self.y_var[start_num+1, t]] + uc_gen_lb[gen_type]*sol[self.y_var[start_num+1, \
                                        t-k]] + sum((uc_gen_lb[gen_type] + (k-t+s)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])* \
                                        (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s-1]]) for s in s_set)
                                    if w > w_bar:
                                        #if the current violation is greater under given t, k, and m, update m_bar
                                        w_bar = w
                                        m_bar = m
                                    m -= 1
                            if w_bar > w_hat:
                                # if current violation under given t and given k is greater than w_hat, update w_hat, k_hat, and m_hat
                                w_hat = w_bar
                                m_hat = m_bar
                                k_hat = k
                        if w_hat > w_most_violate_value:
                            # if w_hat is greater than the current most violation value under given t, update most violation value and indexes
                            w_most_violate_value = w_hat
                            k_most_violate_index = k_hat
                            m_most_violate_index = m_hat
                            t_most_violate_index = t
                            
                    if w_most_violate_value > self.eps:
                        #if the largest violation of (11a) is greater than the violation tolerance epsilon, add it as a constraint to the model
                        alpha = min(k_most_violate_index-1, uc_min_on[gen_type] - 2)
                        s_set = np.array([], dtype = int)
                        for i in range(0, alpha - m_most_violate_index + 1):
                            #find all i in [0, alpha - m*] such that y_[t*-i] - y_[t*-i-1] > 0 to form s_set
                            if sol[self.y_var[start_num+1, t_most_violate_index-i]] - sol[self.y_var[start_num+1, t_most_violate_index-i-1]] > 0:
                                s_set = np.append(s_set, i)
                        for g in range(start_num+1, end_num+1):
                            '''
                            #add the largest violation constraint to all generators of current type
                            self.nb_cuts += 1
                            #update the number of cuts we added
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t_most_violate_index] - self.x_var[g, \
                                    t_most_violate_index - k_most_violate_index] - (uc_gen_lb[gen_type] + (k_most_violate_index - \
                                    m_most_violate_index)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*self.y_var[g, t_most_violate_index \
                                    + m_most_violate_index + 1] - uc_rampup_ub[gen_type]*sum(self.y_var[g, t_most_violate_index + j] \
                                    for j in range(1, m_most_violate_index+1)) - uc_start_ub[gen_type]*self.y_var[g, t_most_violate_index] \
                                    + uc_gen_lb[gen_type]*self.y_var[g, t_most_violate_index - k_most_violate_index] +\
                                    sum((uc_gen_lb[gen_type] + (k_most_violate_index - i)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])* \
                                    (self.y_var[g, t_most_violate_index - i] - self.y_var[g, t_most_violate_index - i - 1]) for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                            '''
                            for t in range(k_most_violate_index+1, uc_opt_period - int(m_most_violate_index)):
                                ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t] - self.x_var[g, \
                                        t - k_most_violate_index] - (uc_gen_lb[gen_type] + (k_most_violate_index - \
                                        m_most_violate_index)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])*self.y_var[g, t\
                                        + m_most_violate_index + 1] - uc_rampup_ub[gen_type]*sum(self.y_var[g, t + j] \
                                        for j in range(1, m_most_violate_index+1)) - uc_start_ub[gen_type]*self.y_var[g, t] \
                                        + uc_gen_lb[gen_type]*self.y_var[g, t - k_most_violate_index] +\
                                        sum((uc_gen_lb[gen_type] + (k_most_violate_index - i)*uc_rampup_ub[gen_type] - uc_start_ub[gen_type])* \
                                        (self.y_var[g, t - i] - self.y_var[g, t - i - 1]) for i in s_set) <= 0)
                                self.add(ct_lhs, ct_sense, ct_rhs)
                    
                #move the start index to the end index of current type of generators
                start_num = end_num
    
    def separation_double_var_rampup_11b(self):
        '''double var exact separation for constraint (11b)'''
        #get the current relaxation solution
        sol = self.make_complete_solution()
        #use start_num and end_num to denote the start index and end index of current type of generators under 
        #current instance
        start_num = 0
        end_num = 0
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #the number of current type of generators for current instance should be greater than 0
                #update the end_num
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                if uc_min_on[gen_type] > 1:
                    # the minimum on period should be greater than 1, i.e., L > 1
                    #define the most violation value, indexes of current relaxation solution
                    w_most_violate_value = 0
                    t_most_violate_index = 0
                    k_most_violate_index = 0
                    m_most_violate_index = 0
                    #calculate the largest violation for each t in [2, T-1]
                    for t in range(2, uc_opt_period):
                        #set the value of k', i.e., the largest possible value of k
                        k_prime = min(uc_opt_period - t, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampdown_ub[gen_type]))
                        #define the largest violation value and index for given t
                        w_hat = 0
                        k_hat = 0
                        m_hat = 0
                        for k in range(1, k_prime+1):
                            #define the upper bound of s_q, i.e., alpha, and the upper bound of m, i.e., beta
                            alpha = min(k-1, uc_min_on[gen_type]-2)
                            beta = min(alpha, t-2)
                            #define the largest violation value and indexes for given t and given k
                            w_bar = 0
                            m_bar = 0
                            s_set = np.array([])
                            #find the value of m0
                            m0 = 0
                            temp0 = -(uc_gen_lb[gen_type] + (k-m0)*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])*sol[self.y_var[start_num+1, \
                                    t-m0-1]] - uc_rampdown_ub[gen_type]*sum(sol[self.y_var[start_num+1, t-j]] for j in range(1, m0+1))
                            for m in range(0, beta+1):
                                temp1 = -(uc_gen_lb[gen_type] + (k-m)*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])* \
                                      sol[self.y_var[start_num+1, t-m-1]] - uc_rampdown_ub[gen_type]*sum(sol[self.y_var[start_num+1, t-j]] for \
                                      j in range(1, m+1))
                                if temp1 > temp0:
                                    temp0 = temp1
                                    m0 = m
                            for i in range(0, alpha-m0+1):
                                #search for i between 0 and alpha - m0 such that y_[t+i] - y_[t+i+1] > 0
                                if sol[self.y_var[start_num+1, t+i]] - sol[self.y_var[start_num+1, t+i+1]] > 0:
                                    s_set = np.append(s_set, t+i)
                            #calculate the violation value for m0
                            w = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t+k]] - (uc_gen_lb[gen_type] + \
                                (k-m0)*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])*sol[self.y_var[start_num+1, t-m0-1]] - \
                                uc_rampdown_ub[gen_type]*sum(sol[self.y_var[start_num+1, t-j]] for j in range(1, m0+1)) - \
                                uc_start_ub[gen_type]*sol[self.y_var[start_num+1, t]] + uc_gen_lb[gen_type]*sol[self.y_var[start_num+1, \
                                t+k]] + sum((uc_gen_lb[gen_type] + (k-s+t)*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])* \
                                (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s+1]]) for s in s_set)
                            if w > w_bar:
                                w_bar = w
                                m_bar = m0
                            if m0 > 0:
                                #if m0 is not 0, we check for any better violation of m in [0,m-1]
                                m = m0-1
                                while m >= 0:
                                    if sol[self.y_var[start_num+1, t+alpha-m]] - sol[self.y_var[start_num+1, t+alpha-m+1]] > 0:
                                        s_set = np.append(s_set, t+alpha-m)
                                    w = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t+k]] - (uc_gen_lb[gen_type] + \
                                        (k-m)*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])*sol[self.y_var[start_num+1, t-m-1]] - \
                                        uc_rampdown_ub[gen_type]*sum(sol[self.y_var[start_num+1, t-j]] for j in range(1, m+1)) - \
                                        uc_start_ub[gen_type]*sol[self.y_var[start_num+1, t]] + uc_gen_lb[gen_type]*sol[self.y_var[start_num+1, \
                                        t+k]] + sum((uc_gen_lb[gen_type] + (k-s+t)*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])* \
                                        (sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s+1]]) for s in s_set)
                                    if w > w_bar:
                                        #if the current violation is greater under given t, k, and m, update m_bar
                                        w_bar = w
                                        m_bar = m
                                    m -= 1
                            if w_bar > w_hat:
                                # if current violation under given t and given k is greater than w_hat, update w_hat, k_hat, and m_hat
                                w_hat = w_bar
                                m_hat = m_bar
                                k_hat = k
                        if w_hat > w_most_violate_value:
                            # if w_hat is greater than the current most violation value under given t, update most violation value and indexes
                            w_most_violate_value = w_hat
                            k_most_violate_index = k_hat
                            m_most_violate_index = m_hat
                            t_most_violate_index = t
                            
                    if w_most_violate_value > self.eps:
                        #if the largest violation of (11a) is greater than the violation tolerance epsilon, add it as a constraint to the model
                        alpha = min(k_most_violate_index-1, uc_min_on[gen_type] - 2)
                        s_set = np.array([], dtype = int)
                        for i in range(0, alpha - m_most_violate_index + 1):
                            #find all i in [0, alpha - m*] such that y_[t*+i] - y_[t*+i+1] > 0 to form s_set
                            if sol[self.y_var[start_num+1, t_most_violate_index+i]] - sol[self.y_var[start_num+1, t_most_violate_index+i+1]] > 0:
                                s_set = np.append(s_set, i)
                        for g in range(start_num+1, end_num+1):
                            '''
                            #add the largest violation constraint to all generators of current type
                            self.nb_cuts += 1
                            #update the number of cuts we added
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t_most_violate_index] - self.x_var[g, \
                                    t_most_violate_index + k_most_violate_index] - (uc_gen_lb[gen_type] + (k_most_violate_index - \
                                    m_most_violate_index)*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])*self.y_var[g, t_most_violate_index \
                                    - m_most_violate_index - 1] - uc_rampdown_ub[gen_type]*sum(self.y_var[g, t_most_violate_index - j] \
                                    for j in range(1, m_most_violate_index+1)) - uc_start_ub[gen_type]*self.y_var[g, t_most_violate_index] \
                                    + uc_gen_lb[gen_type]*self.y_var[g, t_most_violate_index + k_most_violate_index] +\
                                    sum((uc_gen_lb[gen_type] + (k_most_violate_index - i)*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])* \
                                    (self.y_var[g, t_most_violate_index + i] - self.y_var[g, t_most_violate_index + i + 1]) for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                            '''
                            for t in range(int(m_most_violate_index)+2, uc_opt_period-k_most_violate_index+1):
                                ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t] - self.x_var[g, \
                                        t + k_most_violate_index] - (uc_gen_lb[gen_type] + (k_most_violate_index - \
                                        m_most_violate_index)*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])*self.y_var[g, t\
                                        - m_most_violate_index - 1] - uc_rampdown_ub[gen_type]*sum(self.y_var[g, t - j] \
                                        for j in range(1, m_most_violate_index+1)) - uc_start_ub[gen_type]*self.y_var[g, t] \
                                        + uc_gen_lb[gen_type]*self.y_var[g, t + k_most_violate_index] +\
                                        sum((uc_gen_lb[gen_type] + (k_most_violate_index - i)*uc_rampdown_ub[gen_type] - uc_start_ub[gen_type])* \
                                        (self.y_var[g, t + i] - self.y_var[g, t + i + 1]) for i in s_set) <= 0)
                                self.add(ct_lhs, ct_sense, ct_rhs)
                    
                #move the start index to the end index of current type of generators
                start_num = end_num
                
    def separation_double_var_rampup_12a(self):
        '''exact separation for constraint 12a'''
        #get the current relaxation solution
        sol = self.make_complete_solution()
        #use start_num and end_num to denote the start index and end index of current type of generators under 
        #current instance
        start_num = 0
        end_num = 0
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #the number of current type of generators for current instance should be greater than 0
                #update the end_num
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                #define the most violation value, indexes of current relaxation solution
                w_most_violate_value = 0
                t_most_violate_index = 0
                k_most_violate_index = 0
                m_most_violate_index = 0
                #calculate the largest violation for each t in [2, T]
                for t in range(2, uc_opt_period+1):
                    #set the value of k', i.e., the largest possible value of k
                    k_prime = min(t-1, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampup_ub[gen_type]))
                    #define the largest violation value and index for given t
                    w_hat = 0
                    k_hat = 0
                    m_hat = 0
                    for k in range(1, k_prime+1):
                        #define the upper bound of s_q, i.e., alpha, and the upper bound of m, i.e., beta
                        alpha = min(k-1, uc_min_on[gen_type]-1)
                        beta = min(alpha, uc_opt_period - t)
                        #define the largest violation value and indexes for given t and given k
                        w_bar = 0
                        m_bar = 0
                        s_set = np.array([])
                        #find the value of m0
                        m0 = 0
                        temp0 = -uc_rampup_ub[gen_type]*sum(sol[self.y_var[start_num+1, t+j]] for j in range(1, m0+1)) -\
                               (uc_gen_lb[gen_type] + (k-m0)*uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t]]
                        for m in range(0, beta+1):
                            temp1 = -uc_rampup_ub[gen_type]*sum(sol[self.y_var[start_num+1, t+j]] for j in range(1, m+1)) -\
                                   (uc_gen_lb[gen_type] + (k-m)*uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t]]
                            if temp1 > temp0:
                                temp0 = temp1
                                m0 = m
                        for i in range(0, alpha-m0+1):
                            #search for i between 0 and alpha - m0 such that y_[t-i] - y_[t-i-1] > 0
                            if sol[self.y_var[start_num+1, t-i]] - sol[self.y_var[start_num+1, t-i-1]] > 0:
                                s_set = np.append(s_set, t-i)
                        #calculate the violation value for m0
                        w = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t-k]] + temp0 + uc_gen_lb[gen_type]* \
                            sol[self.y_var[start_num+1, t-k]] + sum((uc_gen_lb[gen_type] + (k-t+s)*uc_rampup_ub[gen_type] - \
                            uc_start_ub[gen_type])*(sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s-1]]) \
                            for s in s_set)
                        if w > w_bar:
                            w_bar = w
                            m_bar = m0
                        if m0 > 0:
                            #if m0 is not 0, we check for any better violation of m in [0,m-1]
                            m = m0-1
                            while m >= 0:
                                if sol[self.y_var[start_num+1, t-alpha+m]] - sol[self.y_var[start_num+1, t-alpha+m-1]] > 0:
                                    s_set = np.append(s_set, t-alpha+m)
                                w = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t-k]] - uc_rampup_ub[gen_type]* \
                                    sum(sol[self.y_var[start_num+1, t+j]] for j in range(1, m+1)) - (uc_gen_lb[gen_type] + (k-m)* \
                                    uc_rampup_ub[gen_type])*sol[self.y_var[start_num+1, t]] + uc_gen_lb[gen_type]* \
                                    sol[self.y_var[start_num+1, t-k]] + sum((uc_gen_lb[gen_type] + (k-t+s)*uc_rampup_ub[gen_type]- \
                                    uc_start_ub[gen_type])*(sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s-1]]) \
                                    for s in s_set)
                                if w > w_bar:
                                    #if the current violation is greater under given t, k, and m, update m_bar
                                    w_bar = w
                                    m_bar = m
                                m -= 1
                        if w_bar > w_hat:
                            # if current violation under given t and given k is greater than w_hat, update w_hat, k_hat, and m_hat
                            w_hat = w_bar
                            m_hat = m_bar
                            k_hat = k
                    if w_hat > w_most_violate_value:
                        # if w_hat is greater than the current most violation value under given t, update most violation value and indexes
                        w_most_violate_value = w_hat
                        k_most_violate_index = k_hat
                        m_most_violate_index = m_hat
                        t_most_violate_index = t
                    
                if w_most_violate_value > self.eps:
                    #if the largest violation of (11a) is greater than the violation tolerance epsilon, add it as a constraint to the model
                    alpha = min(k_most_violate_index-1, uc_min_on[gen_type] - 1)
                    s_set = np.array([], dtype = int)
                    for i in range(0, alpha - m_most_violate_index + 1):
                        #find all i in [0, alpha - m*] such that y_[t*-i] - y_[t*-i-1] > 0 to form s_set
                        if sol[self.y_var[start_num+1, t_most_violate_index-i]] - sol[self.y_var[start_num+1, t_most_violate_index-i-1]] > 0:
                            s_set = np.append(s_set, i)
                    for g in range(start_num+1, end_num+1):
                        '''
                        #add the largest violation constraint to all generators of current type
                        self.nb_cuts += 1
                        #update the number of cuts we added
                        ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t_most_violate_index] - self.x_var[g, \
                                t_most_violate_index - k_most_violate_index] - uc_rampup_ub[gen_type]*sum(self.y_var[g, \
                                t_most_violate_index + j] for j in range(1, m_most_violate_index+1)) - (uc_gen_lb[gen_type] + \
                                (k_most_violate_index - m_most_violate_index)*uc_rampup_ub[gen_type])*self.y_var[g, \
                                t_most_violate_index] + uc_gen_lb[gen_type]*self.y_var[g, t_most_violate_index - \
                                k_most_violate_index] + sum((uc_gen_lb[gen_type] + (k_most_violate_index - i)*uc_rampup_ub[gen_type] - \
                                uc_start_ub[gen_type])*(self.y_var[g, t_most_violate_index - i] - self.y_var[g, \
                                t_most_violate_index-i-1]) for i in s_set) <= 0)
                        self.add(ct_lhs, ct_sense, ct_rhs)
                        '''
                        for t in range(int(k_most_violate_index)+1, uc_opt_period-m_most_violate_index+1):
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t] - self.x_var[g, \
                                    t - k_most_violate_index] - uc_rampup_ub[gen_type]*sum(self.y_var[g, \
                                    t + j] for j in range(1, m_most_violate_index+1)) - (uc_gen_lb[gen_type] + \
                                    (k_most_violate_index - m_most_violate_index)*uc_rampup_ub[gen_type])*self.y_var[g, \
                                    t] + uc_gen_lb[gen_type]*self.y_var[g, t - k_most_violate_index] + \
                                    sum((uc_gen_lb[gen_type] + (k_most_violate_index - i)*uc_rampup_ub[gen_type] - \
                                    uc_start_ub[gen_type])*(self.y_var[g, t - i] - self.y_var[g, \
                                    t-i-1]) for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                
                #move the start index to the end index of current type of generators
                start_num = end_num
                
    def separation_double_var_rampup_12b(self):
        '''exact separation for constraint 12b'''
        #get the current relaxation solution
        sol = self.make_complete_solution()
        #use start_num and end_num to denote the start index and end index of current type of generators under 
        #current instance
        start_num = 0
        end_num = 0
        for gen_type in range(0,8):
            if uc_gen_num[self.instance_index-1, gen_type] > 0:
                #the number of current type of generators for current instance should be greater than 0
                #update the end_num
                end_num += uc_gen_num[self.instance_index-1, gen_type]
                #define the most violation value, indexes of current relaxation solution
                w_most_violate_value = 0
                t_most_violate_index = 0
                k_most_violate_index = 0
                m_most_violate_index = 0
                #calculate the largest violation for each t in [1, T-1]
                for t in range(1, uc_opt_period):
                    #set the value of k', i.e., the largest possible value of k
                    k_prime = min(uc_opt_period-t, math.floor((uc_gen_ub[gen_type] - uc_gen_lb[gen_type])/uc_rampdown_ub[gen_type]))
                    #define the largest violation value and index for given t
                    w_hat = 0
                    k_hat = 0
                    m_hat = 0
                    for k in range(1, k_prime+1):
                        #define the upper bound of s_q, i.e., alpha, and the upper bound of m, i.e., beta
                        alpha = min(k-1, uc_min_on[gen_type]-1)
                        beta = min(alpha, t-1)
                        #define the largest violation value and indexes for given t and given k
                        w_bar = 0
                        m_bar = 0
                        s_set = np.array([])
                        #find the value of m0
                        m0 = 0
                        temp0 = -uc_rampdown_ub[gen_type]*sum(sol[self.y_var[start_num+1, t-j]] for j in range(1, m0+1)) -\
                               (uc_gen_lb[gen_type] + (k-m0)*uc_rampdown_ub[gen_type])*sol[self.y_var[start_num+1, t]]
                        for m in range(0, beta+1):
                            temp1 = -uc_rampdown_ub[gen_type]*sum(sol[self.y_var[start_num+1, t-j]] for j in range(1, m+1)) -\
                                   (uc_gen_lb[gen_type] + (k-m)*uc_rampdown_ub[gen_type])*sol[self.y_var[start_num+1, t]]
                            if temp1 > temp0:
                                temp0 = temp1
                                m0 = m
                        for i in range(0, alpha-m0+1):
                            #search for i between 0 and alpha - m0 such that y_[t+i] - y_[t+i+1] > 0
                            if sol[self.y_var[start_num+1, t+i]] - sol[self.y_var[start_num+1, t+i+1]] > 0:
                                s_set = np.append(s_set, t+i)
                        #calculate the violation value for m0
                        w = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t+k]] + temp0 + uc_gen_lb[gen_type]* \
                            sol[self.y_var[start_num+1, t+k]] + sum((uc_gen_lb[gen_type] + (k-s+t)*uc_rampdown_ub[gen_type] - \
                            uc_start_ub[gen_type])*(sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s+1]]) \
                            for s in s_set)
                        if w > w_bar:
                            w_bar = w
                            m_bar = m0
                        if m0 > 0:
                            #if m0 is not 0, we check for any better violation of m in [0,m-1]
                            m = m0-1
                            while m >= 0:
                                if sol[self.y_var[start_num+1, t+alpha-m]] - sol[self.y_var[start_num+1, t+alpha-m+1]] > 0:
                                    s_set = np.append(s_set, t+alpha-m)
                                w = sol[self.x_var[start_num+1, t]] - sol[self.x_var[start_num+1, t+k]] - uc_rampdown_ub[gen_type]* \
                                    sum(sol[self.y_var[start_num+1, t-j]] for j in range(1, m+1)) - (uc_gen_lb[gen_type] + (k-m)* \
                                    uc_rampdown_ub[gen_type])*sol[self.y_var[start_num+1, t]] + uc_gen_lb[gen_type]* \
                                    sol[self.y_var[start_num+1, t+k]] + sum((uc_gen_lb[gen_type] + (k-s+t)*uc_rampdown_ub[gen_type]- \
                                    uc_start_ub[gen_type])*(sol[self.y_var[start_num+1, s]] - sol[self.y_var[start_num+1, s+1]]) \
                                    for s in s_set)
                                if w > w_bar:
                                    #if the current violation is greater under given t, k, and m, update m_bar
                                    w_bar = w
                                    m_bar = m
                                m -= 1
                        if w_bar > w_hat:
                            # if current violation under given t and given k is greater than w_hat, update w_hat, k_hat, and m_hat
                            w_hat = w_bar
                            m_hat = m_bar
                            k_hat = k
                    if w_hat > w_most_violate_value:
                        # if w_hat is greater than the current most violation value under given t, update most violation value and indexes
                        w_most_violate_value = w_hat
                        k_most_violate_index = k_hat
                        m_most_violate_index = m_hat
                        t_most_violate_index = t
                    
                if w_most_violate_value > self.eps:
                    #if the largest violation of (11a) is greater than the violation tolerance epsilon, add it as a constraint to the model
                    alpha = min(k_most_violate_index-1, uc_min_on[gen_type] - 1)
                    s_set = np.array([], dtype = int)
                    for i in range(0, alpha - m_most_violate_index + 1):
                        #find all i in [0, alpha - m*] such that y_[t*-i] - y_[t*-i-1] > 0 to form s_set
                        if sol[self.y_var[start_num+1, t_most_violate_index+i]] - sol[self.y_var[start_num+1, t_most_violate_index+i+1]] > 0:
                            s_set = np.append(s_set, i)
                    for g in range(start_num+1, end_num+1):
                        '''
                        #add the largest violation constraint to all generators of current type
                        self.nb_cuts += 1
                        #update the number of cuts we added
                        ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t_most_violate_index] - self.x_var[g, \
                                t_most_violate_index + k_most_violate_index] - uc_rampdown_ub[gen_type]*sum(self.y_var[g, \
                                t_most_violate_index - j] for j in range(1, m_most_violate_index+1)) - (uc_gen_lb[gen_type] + \
                                (k_most_violate_index - m_most_violate_index)*uc_rampdown_ub[gen_type])*self.y_var[g, \
                                t_most_violate_index] + uc_gen_lb[gen_type]*self.y_var[g, t_most_violate_index + \
                                k_most_violate_index] + sum((uc_gen_lb[gen_type] + (k_most_violate_index - i)*uc_rampdown_ub[gen_type] - \
                                uc_start_ub[gen_type])*(self.y_var[g, t_most_violate_index + i] - self.y_var[g, \
                                t_most_violate_index+i+1]) for i in s_set) <= 0)
                        self.add(ct_lhs, ct_sense, ct_rhs)
                        '''
                        for t in range(int(m_most_violate_index)+1, uc_opt_period-k_most_violate_index+1):
                            ct_lhs, ct_sense, ct_rhs = self.linear_ct_to_cplex(self.x_var[g, t] - self.x_var[g, \
                                    t + k_most_violate_index] - uc_rampdown_ub[gen_type]*sum(self.y_var[g, \
                                    t - j] for j in range(1, m_most_violate_index+1)) - (uc_gen_lb[gen_type] + \
                                    (k_most_violate_index - m_most_violate_index)*uc_rampdown_ub[gen_type])*self.y_var[g, \
                                    t] + uc_gen_lb[gen_type]*self.y_var[g, t + k_most_violate_index] + \
                                    sum((uc_gen_lb[gen_type] + (k_most_violate_index - i)*uc_rampdown_ub[gen_type] - \
                                    uc_start_ub[gen_type])*(self.y_var[g, t + i] - self.y_var[g, \
                                    t+i+1]) for i in s_set) <= 0)
                            self.add(ct_lhs, ct_sense, ct_rhs)
                
                #move the start index to the end index of current type of generators
                start_num = end_num
                            
                            
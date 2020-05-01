# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 11:29:57 2018

@author: jms4
"""
import numpy as np

#%% step through ode
def time_stepper(drive_current_ISI,tau_int,num_input_spikes,index_first_spike,dt,w,voltage_threshold_0,voltage_threshold_delta,tau_ref):   
    time_sim = num_input_spikes*drive_current_ISI+index_first_spike*dt
    num_time_steps = int(time_sim/dt)
    #time_first_spike = index_first_spike*dt
    time_vec = dt*np.linspace(0,num_time_steps-1,num_time_steps)
    print("len(time_vec) = ",len(time_vec))
    num_time_steps_per_input_spike = int(round(drive_current_ISI/dt))
    
    #create vector of spikes
    I_current_drive_vec = np.zeros(num_time_steps)
    for qq in range(num_input_spikes):
        #print(index_first_spike+qq*num_time_steps_per_input_spike-1)
        if np.less_equal(index_first_spike+qq*num_time_steps_per_input_spike-1,len(I_current_drive_vec)-1):
            I_current_drive_vec[index_first_spike+qq*num_time_steps_per_input_spike-1] = w/dt
    
    #solve ODE with Euler method
    voltage_vec = np.zeros(num_time_steps) 
    voltage_threshold_vec = np.zeros(num_time_steps)
    num_spikes = 0
    
    #last spike time
    t_s = -10
    for pp in range(num_time_steps-1):  
        #print("num_spikes = ", num_spikes)
        voltage_vec[pp+1] = (1-dt/tau_int)*voltage_vec[pp]+dt*I_current_drive_vec[pp]
        voltage_threshold = voltage_threshold_0+voltage_threshold_delta*np.exp(-(time_vec[pp]-t_s)/tau_ref)   
        voltage_threshold_vec[pp+1] = voltage_threshold       
        if np.greater_equal(voltage_vec[pp+1],voltage_threshold):
            t_s = time_vec[pp]
            num_spikes += 1
            voltage_vec[pp+1] = 0
    firing_rate = num_spikes/time_sim
    
    print("num_spikes = ", num_spikes)
    print("firing_rate = ", firing_rate*1e-6, " MHz")
    
    return firing_rate,num_spikes
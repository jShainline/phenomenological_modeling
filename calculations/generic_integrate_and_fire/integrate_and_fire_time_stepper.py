# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 11:29:57 2018

@author: jms4
"""
from pylab import *
import numpy as np

#%% step through ode
def time_stepper(drive_current,tau_int,dt,w,voltage_threshold_0,voltage_threshold_delta,tau_ref,time_sim):   
    num_time_steps = int(time_sim/dt)
    time_vec = dt*np.linspace(0,num_time_steps-1,num_time_steps)
        
    #solve ODE with Euler method
    voltage_vec = np.zeros(num_time_steps) 
    voltage_threshold_vec = np.zeros(num_time_steps)
    num_spikes = 0
    
    #last spike time
    t_s = -10
    for pp in range(num_time_steps-1):  
        #print("num_spikes = ", num_spikes)
        voltage_vec[pp+1] = (1-dt/tau_int)*voltage_vec[pp]+dt*w*drive_current
        voltage_threshold = voltage_threshold_0+voltage_threshold_delta*np.exp(-(time_vec[pp]-t_s)/tau_ref)   
        voltage_threshold_vec[pp+1] = voltage_threshold       
        if np.greater_equal(voltage_vec[pp+1],voltage_threshold):
            t_s = time_vec[pp]
            num_spikes += 1
            voltage_vec[pp+1] = 0
    firing_rate = num_spikes/time_sim
    
    print("num_spikes = ", num_spikes)
    print("firing_rate = ", firing_rate*1e-6, " MHz")

#    figure()
#    plot(1e9*time_vec,voltage_vec, label = 'voltage_vec')
#    plot(1e9*time_vec,voltage_threshold_vec, label = 'voltage_threshold_vec')
#    xlabel("time [ns]")
#    legend(loc="best")
#    show()
    
    return firing_rate,num_spikes

def time_stepper_hard_refractory(drive_current,tau_int,dt,w,voltage_threshold_0,tau_ref,time_sim):   
    num_time_steps = int(time_sim/dt)
    time_vec = dt*np.linspace(0,num_time_steps-1,num_time_steps)
        
    #solve ODE with Euler method
    voltage_vec = np.zeros(num_time_steps) 
    num_spikes = 0        
    t_s = -10#last spike time
    for pp in range(num_time_steps-1):  
        #print("num_spikes = ", num_spikes)
        voltage_vec[pp+1] = (1-dt/tau_int)*voltage_vec[pp]+dt*w*drive_current     
        if np.greater_equal(voltage_vec[pp+1],voltage_threshold_0) and np.greater_equal(time_vec[pp],t_s+tau_ref) :
            t_s = time_vec[pp]
            num_spikes += 1
            voltage_vec[pp+1] = 0
    firing_rate = num_spikes/time_sim
    
    print("num_spikes = ", num_spikes)
    print("firing_rate = ", firing_rate*1e-6, " MHz")

#    figure()
#    plot(1e9*time_vec,voltage_vec, label = 'voltage_vec')
#    plot(1e9*time_vec,voltage_threshold_vec, label = 'voltage_threshold_vec')
#    xlabel("time [ns]")
#    legend(loc="best")
#    show()
    
    return firing_rate,num_spikes
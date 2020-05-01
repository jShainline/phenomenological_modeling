#%% import
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 09:58:13 2018

@author: jms4
"""
from pylab import *
#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np
#import math
from integrate_and_fire_time_stepper import time_stepper
from integrate_and_fire_time_stepper import time_stepper_hard_refractory

#%% inputs
numNeurons = 1

drive_current_vec = np.logspace(4,9,300)
#drive_current_vec = np.linspace(1e4,1e9,200)

#time variables
tau_int_vec = np.array([1e-4,1e-5,1e-6,1e-7])#integration time
tau_ref_vec = np.array([50e-9,100e-9])#refractory period

#index_first_spike = 10
#num_input_spikes = 1000
#num_time_steps = 10000
dt = 1e-9

#synaptic weights
w = 1

#threshold
voltage_threshold_vec = np.array([1,10])
#voltage_threshold_delta = 10*voltage_threshold_0

#output data
#num_spikes_mat = np.zeros((num_int_times,num_rates))
#firing_rate_mat = np.zeros((num_int_times,num_rates))

    
firing_rate_mat = np.zeros((len(drive_current_vec),len(tau_int_vec),len(tau_ref_vec),len(voltage_threshold_vec)))
time_sim = dt*100000
for ll in range(len(voltage_threshold_vec)):
    print("ll = ",ll+1," of ",len(voltage_threshold_vec))
    voltage_threshold_0 = voltage_threshold_vec[ll]*w
    for kk in range(len(tau_ref_vec)):
        print("kk = ",kk+1," of ",len(tau_ref_vec))
        tau_ref = tau_ref_vec[kk]
        for jj in range(len(tau_int_vec)):
            print("jj = ",jj+1," of ",len(tau_int_vec))
            tau_int = tau_int_vec[jj]
            for ii in range(len(drive_current_vec)):
                print("ii = ",ii+1," of ",len(drive_current_vec))
                drive_current = drive_current_vec[ii]
                #firing_rate, num_spikes = time_stepper(drive_current,tau_int,dt,w,voltage_threshold_0,voltage_threshold_delta,tau_ref,time_sim)
                firing_rate, num_spikes = time_stepper_hard_refractory(drive_current,tau_int,dt,w,voltage_threshold_0,tau_ref,time_sim)
                firing_rate_mat[ii,jj,kk,ll] = firing_rate

#%% plot

  
#figure()
#plot(drive_current_vec,1e-6*firing_rate_mat[:,0,0], label = 'tau_int = 1e-5,tau_ref = 1e-9')
#plot(drive_current_vec,1e-6*firing_rate_mat[:,0,1], label = 'tau_int = 1e-5,tau_ref = 10e-9')
#plot(drive_current_vec,1e-6*firing_rate_mat[:,0,2], label = 'tau_int = 1e-5,tau_ref = 100e-9')
#plot(drive_current_vec,1e-6*firing_rate_mat[:,1,0], label = 'tau_int = 1e-6,tau_ref = 1e-9')
#plot(drive_current_vec,1e-6*firing_rate_mat[:,1,1], label = 'tau_int = 1e-6,tau_ref = 10e-9')
#plot(drive_current_vec,1e-6*firing_rate_mat[:,1,2], label = 'tau_int = 1e-6,tau_ref = 100e-9')
#plot(drive_current_vec,1e-6*firing_rate_mat[:,2,0], label = 'tau_int = 1e-7,tau_ref = 1e-9')
#plot(drive_current_vec,1e-6*firing_rate_mat[:,2,1], label = 'tau_int = 1e-7,tau_ref = 10e-9')
#plot(drive_current_vec,1e-6*firing_rate_mat[:,2,2], label = 'tau_int = 1e-7,tau_ref = 100e-9')
#xlabel("Drive Current [a.u.]")
#ylabel("Output spike rate [MHz]")
#legend(loc="best")
#show()  
  
figure()
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,0,0,0], label = 't_i=1e-4,t_r=50e-9,v0=1')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,0,0,1], label = 't_i=1e-4,t_r=50e-9,v0=10')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,0,1,0], label = 't_i=1e-4,t_r=100e-9,v0=1')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,0,1,1], label = 't_i=1e-4,t_r=100e-9,v0=10')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,1,0,0], label = 't_i=1e-5,t_r=50e-9,v0=1')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,1,0,1], label = 't_i=1e-5,t_r=50e-9,v0=10')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,1,1,0], label = 't_i=1e-5,t_r=100e-9,v0=1')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,1,1,1], label = 't_i=1e-5,t_r=100e-9,v0=10')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,2,0,0], label = 't_i=1e-6,t_r=50e-9,v0=1')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,2,0,1], label = 't_i=1e-6,t_r=50e-9,v0=10')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,2,1,0], label = 't_i=1e-6,t_r=100e-9,v0=1')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,2,1,1], label = 't_i=1e-6,t_r=100e-9,v0=10')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,3,0,0], label = 't_i=1e-7,t_r=50e-9,v0=1')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,3,0,1], label = 't_i=1e-7,t_r=50e-9,v0=10')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,3,1,0], label = 't_i=1e-7,t_r=100e-9,v0=1')
semilogx(drive_current_vec,1e-6*firing_rate_mat[:,3,1,1], label = 't_i=1e-7,t_r=100e-9,v0=10')
xlabel("Drive Current [a.u.]")
ylabel("Output spike rate [MHz]")
legend(loc="best")
show() 
  

#figure()
#loglog(drive_current_vec,1e-6*firing_rate_mat[:,0], label = 'tau_int_vec[0]')
#loglog(drive_current_vec,1e-6*firing_rate_mat[:,1], label = 'tau_int_vec[1]')
#loglog(drive_current_vec,1e-6*firing_rate_mat[:,2], label = 'tau_int_vec[2]')
#xlabel("Drive Current [a.u.]")
#ylabel("Output spike rate [MHz]")
#legend(loc="best")
#show() 

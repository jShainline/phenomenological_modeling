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

#%% inputs
numNeurons = 1

#drive_current_ISI_vec = np.logspace(-9,-6,100)
drive_current_ISI_vec = np.linspace(1e-8,1e-6,20)

#time variables
tau_int_vec = np.array([1e-5])#np.array([1e-5,1e-6,1e-7])#integration time
tau_ref = 50e-9#refractory period

index_first_spike = 10
num_input_spikes = 100
#num_time_steps = 10000
dt = 1e-9

#synaptic weights
w = 10

#threshold
voltage_threshold_0 = 10*w
voltage_threshold_delta = 10*voltage_threshold_0


#output data
#num_spikes_mat = np.zeros((num_int_times,num_rates))
#firing_rate_mat = np.zeros((num_int_times,num_rates))


firing_rate_mat = np.zeros((len(drive_current_ISI_vec),len(tau_int_vec)))
for jj in range(len(tau_int_vec)):
    print("jj = ",jj+1," of ",len(tau_int_vec))
    tau_int = tau_int_vec[jj]
    for ii in range(len(drive_current_ISI_vec)):
        print("ii = ",ii+1," of ",len(drive_current_ISI_vec))
        drive_current_ISI = drive_current_ISI_vec[ii]
        firing_rate, num_spikes = time_stepper(drive_current_ISI,tau_int,num_input_spikes,index_first_spike,dt,w,voltage_threshold_0,voltage_threshold_delta,tau_ref)
        firing_rate_mat[ii,jj] = firing_rate

#%% plot
drive_current_rate = np.flip(1./drive_current_ISI_vec,0)
  
figure()
plot(drive_current_ISI_vec,firing_rate_mat[:,0], label = 'tau_int_vec[0]')
xlabel("Drive Current ISI [s]")
legend(loc="best")
show()  

temp_vec = firing_rate_mat[:,0]
figure()
plot(drive_current_rate,temp_vec, label = 'tau_int_vec[0]')
xlabel("Drive Current Spike Rate [Hz]")
legend(loc="best")
show()
 
temp_vec = np.flip(firing_rate_mat[:,0],0)
figure()
plot(drive_current_rate,temp_vec, label = 'tau_int_vec[0]')
xlabel("Drive Current Spike Rate [Hz]")
legend(loc="best")
show()

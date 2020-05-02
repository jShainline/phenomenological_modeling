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

#%% inputs
numNeurons = 1

drive_current_ISI = 1e-6

#time variables
tau_int = 1e-5#integration time
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

#last spike time
t_s = -10

#output data
#num_spikes_mat = np.zeros((num_int_times,num_rates))
#firing_rate_mat = np.zeros((num_int_times,num_rates))

#%% step through ode   
time_sim = num_input_spikes*drive_current_ISI+index_first_spike*dt
num_time_steps = int(time_sim/dt)
time_first_spike = index_first_spike*dt
time_vec = [dt*qq for qq in range(num_time_steps)]
num_time_steps_per_input_spike = int(round(drive_current_ISI/dt))

#create vector of spikes
I_current_drive_vec = [0]*num_time_steps
for qq in range(num_input_spikes):
    #print(index_first_spike+qq*num_time_steps_per_input_spike-1)
    I_current_drive_vec[index_first_spike+qq*num_time_steps_per_input_spike-1] = w/dt

#solve ODE with Euler method
voltage_vec = [0]*num_time_steps 
voltage_threshold_vec = [0]*num_time_steps
num_spikes = 0
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

#%% plot
figure()
time_vec = [x*1e9 for x in time_vec]
plot(time_vec,I_current_drive_vec, label = 'Current Drive')
plot(time_vec,voltage_vec, label = 'Membrane Voltage')
xlabel("time [ns]")
legend(loc="best")
show()

figure()
plot(time_vec,I_current_drive_vec, label = 'Current Drive')
xlabel("time [ns]")
legend(loc="best")
show()

figure()
plot(time_vec,voltage_vec, label = 'Membrane Voltage')
xlabel("time [ns]")
legend(loc="best")
show()
#drive_current_rate_vec = 1./drive_current_ISI_vec
#ax.plot(drive_current_rate_vec,firing_rate_mat[0,:])

#ax.set(xlabel='Rate of input spikes (Hz)', ylabel='Rate of output spikes (Hz)',
       #title='Neuron rate transfer function')
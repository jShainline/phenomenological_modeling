function [I_b] = f_synapse_drive_def(t,input_spike_times,I_b_0,I_0,tau_spd)

I_b = zeros(size(t));

for jj = 1:length(input_spike_times)
    ind_vec = find( t >= input_spike_times(jj) );
    I_b(ind_vec) = I_0*exp( -( t(ind_vec)-input_spike_times(jj) )/tau_spd );
end

I_b = I_b+I_b_0;

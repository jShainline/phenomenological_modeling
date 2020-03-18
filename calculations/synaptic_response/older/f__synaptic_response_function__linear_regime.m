function [I_si_vec] = f__synaptic_response_function__linear_regime(time_vec,input_spike_times,I_0,tau_rise,tau_fall)

I_si_mat = zeros(length(time_vec),length(input_spike_times));
for ii = 1:length(input_spike_times)
    ind_vec = find( time_vec > input_spike_times(ii) );
    I_si_mat(ind_vec(1:end),ii) = I_0*(1-exp(-(time_vec(ind_vec(1:end))-input_spike_times(ii))/tau_rise)).*exp(-(time_vec(ind_vec(1:end))-input_spike_times(ii))/tau_fall);
end
I_si_vec = sum(I_si_mat,2);
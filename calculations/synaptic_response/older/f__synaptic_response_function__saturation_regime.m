function [I_si_vec] = f__synaptic_response_function__saturation_regime(time_vec,input_spike_times,I_0,tau_rise,tau_fall)

I_si_vec = zeros(size(time_vec));
for ii = 1:length(input_spike_times)
    ind_vec = find( time_vec > input_spike_times(ii) );
    I_si_vec(ind_vec(3:end)) = I_0*(1-exp(-(time_vec(ind_vec(3:end))-input_spike_times(ii))/tau_rise)).*exp(-(time_vec(ind_vec(3:end))-input_spike_times(ii))/tau_fall);
end
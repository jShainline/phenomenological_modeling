function [I_si_vec] = f__synaptic_response_function(time_vec,input_spike_times,I_0,I_si_sat,gamma1,gamma2,gamma3,tau_rise,tau_fall)

I_si_mat = zeros(length(time_vec),length(input_spike_times));

for ii = 1:length(input_spike_times)
    
    ind_vec = find( time_vec > input_spike_times(ii) );
    I_si_vec_temp = sum(I_si_mat,2);
    
%         I_si_mat(ind_vec(1:end),ii) = f__synaptic_response_prefactor(I_0,I_si_sat,gamma1,gamma2,I_si_vec_temp(ind_vec(ii)-1))*...
%             (1-exp(-(time_vec(ind_vec(1:end))-input_spike_times(ii))/tau_rise)).*exp(-(time_vec(ind_vec(1:end))-input_spike_times(ii))/tau_fall);
    
    
    for jj = 1:length(ind_vec)
        
        if time_vec(ind_vec(jj))-input_spike_times(ii) <= tau_rise
            
            I_si_mat(ind_vec(jj),ii) = f__synaptic_response_prefactor(I_0,I_si_sat,gamma1,gamma2,I_si_vec_temp(ind_vec(jj)-1),tau_rise,tau_fall)*...
                ( (1/tau_rise^gamma3).*( time_vec(ind_vec(jj)) - input_spike_times(ii) )^gamma3 )*...
                exp(tau_rise/tau_fall).*...
                exp(-(time_vec(ind_vec(jj))-input_spike_times(ii))/tau_fall);
            
        else
            
            I_si_mat(ind_vec(jj),ii) = f__synaptic_response_prefactor(I_0,I_si_sat,gamma1,gamma2,I_si_vec_temp(ind_vec(jj)-1),tau_rise,tau_fall)*...
                exp(tau_rise/tau_fall).*...
                exp(-(time_vec(ind_vec(jj))-input_spike_times(ii))/tau_fall);
            
        end
        
    end
    
    I_si_vec = sum(I_si_mat,2);
    
end

end

function [I_prefactor] = f__synaptic_response_prefactor(I_0,I_si_sat,gamma1,gamma2,I_si,tau_rise,tau_fall)

if I_si >= 0 && I_si < I_si_sat %
    A = I_0;
    I_prefactor = min([A*(1-(I_si/I_si_sat)^gamma1)^gamma2 (I_si_sat-I_si)*exp(tau_rise/tau_fall)]);
    %     I_prefactor = A*(1-log(I_si/I_si_sat)/log(gamma1))^gamma2;
    
    %     fprintf('\n\nI_si = %g',I_si)
    %     fprintf('\n\nI_prefactor = %g',I_prefactor)
    
    %     I_prefactor = I_0*(I_si_sat-I_si)/I_si_sat;
    %     I_prefactor = I_0*(1-exp((I_si_sat-I_si)/I_si_sat));
else
    I_prefactor = 0;
end

end
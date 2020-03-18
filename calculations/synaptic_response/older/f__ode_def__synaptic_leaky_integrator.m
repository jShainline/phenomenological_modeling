function [dI_si_dt] = f__ode_def__synaptic_leaky_integrator(t,I_si,L_si,r_si,Ic_jj,r_jj,L_jj,input_spike_times,I_b_0,I_0,tau_spd)

[I_b] = f_synapse_drive_def(t,input_spike_times,I_b_0,I_0,tau_spd);
dIb_dt = [0 diff(I_b)./diff(t)];

% if ( (L_si/(L_si+L_jj))*I_b-I_si )^2 > Ic_jj^2
%     temp_num = (r_jj/L_si)*( ( (L_si/(L_si+L_jj))*I_b-I_si )^2 - Ic_jj^2 )^(1/2);
% else
%     temp_num = 0;
% end
% dI_si_dt(1) = temp_num - (r_si/L_si)*I_si + (L_jj/(L_si+L_jj))*dIb_dt;

if ( (L_si/(L_si+L_jj))*I_b-I_si )^2 > Ic_jj^2
    temp_num = (r_jj/L_si)*( ( (L_si/(L_si+L_jj))*I_b-I_si )^2 - Ic_jj^2 )^(1/2);
else
    temp_num = 0;
end
dI_si_dt(1) = temp_num - (r_si/L_si)*I_si + (L_jj/(L_si+L_jj))*dIb_dt;
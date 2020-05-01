function [I_tTron,I_tTron_e,I_tTron_i,spikeTimes] = f_odeDef_integrateAndFireModel_refractoryRise_withInhibition(I_tTron_0,ItTronCrit,Imax_i,alpha,Ispd_e,Ispd_i,tau_integrate_e,tau_integrate_i,tau_refractory,jPhDrive_e,jPhDrive_i,tVec,refractoryAffectsInhibition)

nT = length(tVec);
dT = tVec(2)-tVec(1);

I_tTron = zeros(nT,1);
I_tTron_i = I_tTron;
I_tTron_e = I_tTron;

spikeTimes = -1000;

I_tTron(1) = I_tTron_0;
I_tTron_e(1) = 0;
I_tTron_i(1) = Imax_i;


if ~strcmp(refractoryAffectsInhibition,'yes') && ~strcmp(refractoryAffectsInhibition,'no'); error('refractoryAffectsInhibition has to be either yes or no'); end

for ii = 2:nT
    
    % without refractory reset on inhibitory array
    if strcmp(refractoryAffectsInhibition,'no')
        I_tTron_i(ii) = max([0 I_tTron_i(ii-1)+(-alpha*jPhDrive_i(ii)*Ispd_i+((Imax_i-I_tTron_i(ii-1))/tau_integrate_i) )*dT]);
    end
    
    if I_tTron(ii-1) < ItTronCrit
        
        I_tTron_e(ii) = I_tTron_e(ii-1) + ( alpha*jPhDrive_e(ii)*Ispd_e*(1-exp(-(tVec(ii)-spikeTimes(end))/tau_refractory))-(I_tTron_e(ii-1)/tau_integrate_e) )*dT;
        
        % with refractory reset on inhibitory array
        if strcmp(refractoryAffectsInhibition,'yes')
            I_tTron_i(ii) = max([0 I_tTron_i(ii-1)+(-alpha*jPhDrive_i(ii)*Ispd_i*(1-exp(-(tVec(ii)-spikeTimes(end))/tau_refractory))+((Imax_i-I_tTron_i(ii-1))/tau_integrate_i) )*dT]);
        end
        
        I_tTron(ii) = I_tTron_e(ii) + I_tTron_i(ii);
        
    elseif I_tTron(ii-1) >= ItTronCrit
        spikeTimes(end+1) = tVec(ii);
        %         I_tTron_e(ii) = 0;
        %         I_tTron_i(ii) = 0;
    end
end





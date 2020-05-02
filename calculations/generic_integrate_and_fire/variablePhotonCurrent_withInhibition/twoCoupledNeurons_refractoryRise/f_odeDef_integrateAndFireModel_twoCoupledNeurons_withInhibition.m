function [I_tTron,I_tTron_e,I_tTron_i,spikeTimes,jPhDrive_e,jPhCouple_e,jPhCouple_i] = f_odeDef_integrateAndFireModel_twoCoupledNeurons_withInhibition(I_tTron_0,ItTronCrit,Imax_i,alpha,Ispd_e,Ispd_i,tau_integrate_e,tau_integrate_i,tau_refractory,jPhDrive_e,jPhDrive_i,tVec,refractoryAffectsInhibition,photonsPerSpike,typeOfCoupling,spikeTimeDelay,spikeDuration)

nT = length(tVec);
dT = tVec(2)-tVec(1);
spikeDuration_numSteps = round(spikeDuration/dT);
spikeTimeDelay_numSteps = round(spikeTimeDelay/dT);

% initialize
I_tTron = cell(2); I_tTron_i = cell(2); I_tTron_e = cell(2); jPhCouple_e = cell(2); jPhCouple_i = cell(2);
I_tTron{1} = zeros(nT,1); I_tTron{2} = zeros(nT,1);
I_tTron_i{1} = I_tTron{1}; I_tTron_i{2} = I_tTron{1}; I_tTron_e{1} = I_tTron{1}; I_tTron_e{2} = I_tTron{1};
jPhCouple_e{1} = I_tTron{1}; jPhCouple_e{2} = I_tTron{1}; jPhCouple_i{1} = I_tTron{1}; jPhCouple_i{2} = I_tTron{1};


spikeTimes = cell(2);
spikeTimes{1} = -1000; spikeTimes{2} = -1000;

% initial conditions
I_tTron{1}(1) = I_tTron_0; I_tTron{2}(1) = I_tTron_0;
I_tTron_e{1}(1) = 0; I_tTron_e{2}(1) = 0;
I_tTron_i{1}(1) = Imax_i; I_tTron_i{2}(1) = Imax_i;

if ~strcmp(refractoryAffectsInhibition,'yes') && ~strcmp(refractoryAffectsInhibition,'no'); error('refractoryAffectsInhibition has to be either yes or no'); end
if ~strcmp(typeOfCoupling,'excitatory') && ~strcmp(typeOfCoupling,'inhibitory') && ~strcmp(typeOfCoupling,'none'); error('typeOfCoupling has to be excitatory, inhibitory, or none'); end

for ii = 2:nT-1
    
    % without refractory reset on inhibitory array
    if strcmp(refractoryAffectsInhibition,'no')
        
        I_tTron_i{1}(ii) = max([0 I_tTron_i{1}(ii-1)+( -alpha*(jPhDrive_i{1}(ii)+jPhCouple_i{1}(ii))*Ispd_i+((Imax_i-I_tTron_i{1}(ii-1))/tau_integrate_i) )*dT]);
        I_tTron_i{2}(ii) = max([0 I_tTron_i{2}(ii-1)+( -alpha*(jPhDrive_i{2}(ii)+jPhCouple_i{2}(ii))*Ispd_i+((Imax_i-I_tTron_i{2}(ii-1))/tau_integrate_i) )*dT]);
        
    end
    
    if I_tTron{1}(ii-1) < ItTronCrit
        
        I_tTron_e{1}(ii) = I_tTron_e{1}(ii-1) + ( alpha*(jPhDrive_e{1}(ii)+jPhCouple_e{1}(ii))*Ispd_e*(1-exp(-(tVec(ii)-spikeTimes{1}(end))/tau_refractory))-(I_tTron_e{1}(ii-1)/tau_integrate_e) )*dT;
        
        % with refractory reset on inhibitory array
        if strcmp(refractoryAffectsInhibition,'yes')
            I_tTron_i{1}(ii) = max([0 I_tTron_i{1}(ii-1)+(-alpha*(jPhDrive_i{1}(ii)+jPhCouple_i{1}(ii))*Ispd_i*(1-exp(-(tVec(ii)-spikeTimes(end))/tau_refractory))+((Imax_i-I_tTron_i{1}(ii-1))/tau_integrate_i))*dT]);
        end
        
        I_tTron{1}(ii) = I_tTron_e{1}(ii) + I_tTron_i{1}(ii);
        
    elseif I_tTron{1}(ii-1) >= ItTronCrit
        
        spikeTimes{1}(end+1) = tVec(ii);
        
        index1 = ii+spikeTimeDelay_numSteps+1;
        index2 = min([nT ii+spikeTimeDelay_numSteps+spikeDuration_numSteps+1]);
        if strcmp(typeOfCoupling,'excitatory')
            jPhCouple_e{2}(index1:index2) = jPhCouple_e{2}(index1:index2)+photonsPerSpike/spikeDuration_numSteps/dT;
        elseif strcmp(typeOfCoupling,'inhibitory')
            jPhCouple_i{2}(index1:index2) = jPhCouple_i{2}(index1:index2)+photonsPerSpike/spikeDuration_numSteps/dT;
        end
        
    end
    
    if I_tTron{2}(ii-1) < ItTronCrit
        
        I_tTron_e{2}(ii) = I_tTron_e{2}(ii-1) + ( alpha*(jPhDrive_e{2}(ii)+jPhCouple_e{2}(ii))*Ispd_e*(1-exp(-(tVec(ii)-spikeTimes{2}(end))/tau_refractory))-(I_tTron_e{2}(ii-1)/tau_integrate_e) )*dT;
        
        % with refractory reset on inhibitory array
        if strcmp(refractoryAffectsInhibition,'yes')
            I_tTron_i{2}(ii) = max([0 I_tTron_i{2}(ii-1)+(-alpha*(jPhDrive_i{2}(ii)+jPhCouple_i{2}(ii))*Ispd_i*(1-exp(-(tVec(ii)-spikeTimes(end))/tau_refractory))+((Imax_i-I_tTron_i{2}(ii-1))/tau_integrate_i))*dT]);
        end
        
        I_tTron{2}(ii) = I_tTron_e{2}(ii) + I_tTron_i{2}(ii);
        
    elseif I_tTron{2}(ii-1) >= ItTronCrit
        
        spikeTimes{2}(end+1) = tVec(ii);
        
        index1 = ii+spikeTimeDelay_numSteps+1;
        index2 = min([nT ii+spikeTimeDelay_numSteps+spikeDuration_numSteps+1]);
        if strcmp(typeOfCoupling,'excitatory')
            jPhCouple_e{1}(index1:index2) = jPhCouple_e{1}(index1:index2)+photonsPerSpike/spikeDuration_numSteps/dT;
        elseif strcmp(typeOfCoupling,'inhibitory')
            jPhCouple_i{1}(index1:index2) = jPhCouple_i{1}(index1:index2)+photonsPerSpike/spikeDuration_numSteps/dT;
        end
        
    end
    
end





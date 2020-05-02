function [dInTron,InTron,spikeTimes,jPhDrive] = f_odeDef_integrateAndFireModel(InTron_0,dInTron_0,InTronCrit,alpha,Ispd,tau_integrate,tau_refractory,jPhDrive,tVec)

nT = length(tVec);
dT = tVec(2)-tVec(1);

InTron = zeros(nT,1);
dInTron = InTron;
spikeTimes = [];

InTron(1) = InTron_0;
dInTron(1) = dInTron_0;

refractoryCounter = 0;
for ii = 2:nT
    if refractoryCounter == 1;
        if tVec(ii) > spikeTimes(end)+tau_refractory
            refractoryCounter = 0;
        end
    elseif InTron(ii-1) < InTronCrit && refractoryCounter == 0
        dInTron(ii) = (alpha*jPhDrive(ii)*Ispd - InTron(ii-1)/tau_integrate)*dT;
        InTron(ii) = InTron(ii-1)+dInTron(ii);
    elseif InTron(ii-1) >= InTronCrit && refractoryCounter == 0        
        spikeTimes(end+1) = tVec(ii);
        InTron(ii) = 0;
        refractoryCounter = 1;        
    end
end





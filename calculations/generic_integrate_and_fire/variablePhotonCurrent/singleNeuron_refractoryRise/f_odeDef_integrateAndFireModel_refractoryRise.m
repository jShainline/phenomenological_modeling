function [InTron,spikeTimes,jPhDrive] = f_odeDef_integrateAndFireModel_refractoryRise(InTron_0,InTronCrit,alpha,Ispd,tau_integrate,tau_refractory,jPhDrive,tVec)

nT = length(tVec);
dT = tVec(2)-tVec(1);

InTron = zeros(nT,1);
spikeTimes = -1000;

InTron(1) = InTron_0;

for ii = 2:nT
    if InTron(ii-1) < InTronCrit        
%         InTron(ii) = InTron(ii-1)+( alpha*jPhDrive(ii)*Ispd/tau_refractory - InTron(ii-1)/tau_integrate )*dT;
        InTron(ii) = InTron(ii-1)+( alpha*jPhDrive(ii)*Ispd* ( 1 - exp(-(tVec(ii)-spikeTimes(end))/tau_refractory) ) - InTron(ii-1)/tau_integrate )*dT;
    elseif InTron(ii-1) >= InTronCrit        
        spikeTimes(end+1) = tVec(ii);
        InTron(ii) = 0;       
    end
end





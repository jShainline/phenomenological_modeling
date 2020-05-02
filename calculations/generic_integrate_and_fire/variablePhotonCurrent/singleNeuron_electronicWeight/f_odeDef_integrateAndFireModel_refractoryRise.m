function [IyTron,spikeTimes,jPhDrive] = f_odeDef_integrateAndFireModel_refractoryRise(IyTron_0,IyTronCrit,alpha,Ispd,tau_integrate,tau_refractory,jPhDrive,tVec)

nT = length(tVec);
dT = tVec(2)-tVec(1);

IyTron = zeros(nT,1);
spikeTimes = -1000;

IyTron(1) = IyTron_0;

for ii = 2:nT
    if IyTron(ii-1) < IyTronCrit        
%         InTron(ii) = InTron(ii-1)+( alpha*jPhDrive(ii)*Ispd/tau_refractory - InTron(ii-1)/tau_integrate )*dT;
        IyTron(ii) = IyTron(ii-1)+( alpha*jPhDrive(ii)*Ispd* ( 1 - exp(-(tVec(ii)-spikeTimes(end))/tau_refractory) ) - IyTron(ii-1)/tau_integrate )*dT;
    elseif IyTron(ii-1) >= IyTronCrit        
        spikeTimes(end+1) = tVec(ii);
        IyTron(ii) = 0;       
    end
end





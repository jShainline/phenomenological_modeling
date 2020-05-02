function [dInTron,InTron,spikeTimes,jPhDrive_weighted] = f_odeDef_integrateAndFireModel_threeByThreeClassifier(InTron_0,dInTron_0,InTronCrit,alpha,Ispd,tau_integrate,tau_refractory,jPh0Vec,W,tVec)

nT = length(tVec);
dT = tVec(2)-tVec(1);

InTron = zeros(3,nT);%one for each output neuron
dInTron = InTron;
spikeTimes = cell(3,1);%one for each output neuron
for ii = 1:3
    spikeTimes{ii} = -1000;
end

InTron(1) = InTron_0;
dInTron(1) = dInTron_0;

% calculate inputs to each output neuron
% fprintf('\n\nsize(W) = %g\n\n',size(W))
% fprintf('\n\nsize(jPh0Vec) = %g\n\n',size(jPh0Vec))
jPhDrive_weighted = W*jPh0Vec';
% fprintf('\n\nsize(jPhDrive_weighted) = %g\n\n',size(jPhDrive_weighted))

for jj = 1:3%loop over output neurons
    
    for ii = 2:nT
        if InTron(jj,ii-1) < InTronCrit
            dInTron(jj,ii) = (alpha*jPhDrive_weighted(jj)*Ispd* ( 1 - exp(-(tVec(ii)-spikeTimes{jj}(end))/tau_refractory) ) - InTron(jj,ii-1)/tau_integrate)*dT;
            InTron(jj,ii) = InTron(jj,ii-1)+dInTron(jj,ii);
        elseif InTron(jj,ii-1) >= InTronCrit
            spikeTimes{jj}(end+1) = tVec(ii);
            InTron(jj,ii) = 0;
        end
    end
    
end





%% initialize
clc
clear all
close all

[fontName,fontSize,bRGY,scrsz] = f_plotting;

p = f_physicalConstants;

%% input parameters

%parameters
Ispd = 1e-6;%current through a single snspd of the receiver array
jPh0 = 200e6;%photons per second when pulse is on
alpha = 1;%efficiency factor scaling sent photons to received photons
tau_integrate = 1e-6;%receiver integration time
tau_refractory = 50e-9;%refractory period
InTronCrit = Ispd*10;

%temporal
t0 = 0;%beginning of simulation
tf = 1e-6;%total simulation time in seconds
nT = 1e3;%number of time steps
tVec = linspace(t0,tf,nT);

%input signal
tau_switch = 1e-6;
switchTimes = [0.2*tau_switch 1.2*tau_switch];%[timeOn timeOff]; not the same use as some of these other codes; this is ad-hoc to a single square pulse

%inital conditions
InTron_0 = 0;
dInTron_0 = 0;

%% calculate rate response function and its derivative
jPhVec = linspace(1e5,1e10,1e4);
tISI = -tau_integrate * log( 1 - (InTronCrit) ./ (tau_integrate*alpha*jPhVec*Ispd) ) + tau_refractory;
ind = find(real(tISI) > 4e-6);
tISI(ind) = 1;
rISI = 1./tISI;
drISI = diff(rISI);
djPhVec = diff(jPhVec);
drISI_djPh = drISI./djPhVec;

% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(jPhVec,tISI,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
% ylabel('tISI [s]','FontSize',fontSize,'FontName','Times')
% xlabel('jPh [photons/s]','FontSize',fontSize,'FontName','Times')
% % title(sprintf('InTron/InTronCrit for output neuron %g',jj),'FontSize',fontSize,'FontName','Times')
% set(gca,'FontSize',fontSize,'FontName',fontName)
% ylim([0 4e-6])

% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(jPhVec,rISI,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
% ylabel('rISI [Hz]','FontSize',fontSize,'FontName','Times')
% xlabel('jPh [photons/s]','FontSize',fontSize,'FontName','Times')
% % title(sprintf('InTron/InTronCrit for output neuron %g',jj),'FontSize',fontSize,'FontName','Times')
% set(gca,'FontSize',fontSize,'FontName',fontName)

%% training
% initialize weight matrix
[W] = f_createWeightMatrix;

letterList = {'z','v','n'};%'z', 'v', or 'n'
defectVec = 0:9;%0-9; 0 means no defect; 1-9 labels pixel that gets switched; indexing goes down columns

outputRateBook = zeros(3,length(letterList),length(defectVec));
jPhInputBook = zeros(9,length(letterList),length(defectVec));
jPhDriveBook = outputRateBook;
maxRateMat = zeros(length(letterList),length(defectVec));

% letter = 'z';
% defect = 0;
% [jPh0Vec] = f_createInputImageDrive(letter,defect,jPh0);
% [dInTron_dt,InTron,spikeTimes,jPhDrive_weighted] = f_odeDef_integrateAndFireModel_threeByThreeClassifier(InTron_0,dInTron_0,InTronCrit,alpha,Ispd,tau_integrate,tau_refractory,jPh0Vec,W,tVec);
% return

%calcualate maxRate to determine target value
for ii = 1:length(letterList)
    for jj = 1:length(defectVec)
        [jPh0Vec] = f_createInputImageDrive(letterList{ii},defectVec(jj),jPh0);
%         maxInput = sum(jPh0Vec);
maxInput = jPh0;
        [~,ind] = min(abs(jPhVec-maxInput));
        maxRateMat(ii,jj) = rISI(ind);
    end
end

fTarget = zeros(3,length(letterList),length(defectVec));
for ii = 1:3
    for jj = 1:length(letterList)
        if ii == jj
            for kk = 1:length(defectVec)
                fTarget(ii,jj,kk) = maxRateMat(jj,kk);
            end
        end
    end
end

numEpochs = 100;
eta = 0.1;%multiplier that determines rate of learning
for kk = 1:numEpochs
    
    fprintf('\n\nEpoch number %g of %g\n\n',kk,numEpochs)
    
    W
    
    for ii = 1:length(letterList)
        for jj = 1:length(defectVec)
            
            letter = letterList{ii};
            defect = defectVec(jj);
            
            [jPh0Vec] = f_createInputImageDrive(letter,defect,jPh0);
            jPhInputBook(:,ii,jj) = jPh0Vec;
            %             [jPhDrive] = f_photonDriveDef_threeByThreeClassifier(tVec,switchTimes,jPh0Vec);
            % figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
            % plot(tVec*1e6,jPhDrive(1,:),'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
            % hold on
            % plot(tVec*1e6,jPhDrive(2,:),'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
            % ylabel('jPhDrive [photons/sec]','FontSize',fontSize,'FontName','Times')
            % xlabel('Time [\mu s]','FontSize',fontSize,'FontName','Times')
            % legend('pixel 1','pixel 2')
            % % title(sprintf('Total energy: %g fJ\nLED energy: %g fJ\nfraction = %g\ntotalPhotons = %g, assuming efficiency = %g',totEnergy*1e15,ledEnergy*1e15,ledEnergy/totEnergy,totalPhotons,eta),'FontSize',fontSize,'FontName','Times')
            % set(gca,'FontSize',fontSize,'FontName',fontName)
            [dInTron_dt,InTron,spikeTimes,jPhDrive_weighted] = f_odeDef_integrateAndFireModel_threeByThreeClassifier(InTron_0,dInTron_0,InTronCrit,alpha,Ispd,tau_integrate,tau_refractory,jPh0Vec,W,tVec);
            for mm = 1:3
                jPhDriveBook(mm,ii,jj) = jPhDrive_weighted(mm);
                if length(spikeTimes{mm}(:)) > 1
                    numSpikes = length(spikeTimes{mm}(2:end));%throw away the first dummy spike in the distant past
                    diffSpikeTimes = diff(spikeTimes{mm}(2:end));%throw away the first dummy spike in the distant past                    
                    averageInterSpikeInterval = sum(diffSpikeTimes)/numSpikes;
                    averageSpikeRate = 1/averageInterSpikeInterval;
                    outputRateBook(mm,ii,jj) = averageSpikeRate;
                else
                    outputRateBook(mm,ii,jj) = 0;
                end
            end
            
        end
    end
    
    % update weight matrix
    delta = zeros(3,length(letterList),length(defectVec));
    Delta = zeros(3,9,length(letterList),length(defectVec));
    for pp = 1:3%loop over outputs
        for qq = 1:9%loop over inputs
            for rr = 1:length(letterList)%loop over letters
                for ss = 1:length(defectVec)%loop over defects
                    [~,ind] = min( abs( jPhVec-jPhDriveBook(pp,rr,ss) ));
                    dfdI = drISI_djPh(ind);
                    prefactor = fTarget(pp,rr,ss)-outputRateBook(pp,rr,ss);
                    delta(pp,rr,ss) = prefactor*dfdI;
                end
            end
            Delta(pp,qq,rr,ss) = delta(pp,rr,ss)*jPhInputBook(qq,rr,ss);
        end
    end
    
    for pp = 1:3%loop over outputs
        for qq = 1:9%loop over inputs
            summation = 0;
            for rr = 1:length(letterList)%loop over letters
                for ss = 1:length(defectVec)%loop over defects
                    summation = summation+Delta(pp,qq,rr,ss);
                end
            end
            W(pp,qq) = W(pp,qq)+eta*sign(summation);
        end
    end
    
end

W

%% test
correctAnswers = [ones(1,10) 2*ones(1,10) 3*ones(1,10)];
for ii = 1:length(letterList)
    for jj = 1:length(defectVec)
        
        letter = letterList{ii};
        defect = defectVec(jj);
        
        [jPh0Vec] = f_createInputImageDrive(letter,defect,jPh0);
        [dInTron_dt,InTron,spikeTimes,jPhDrive_weighted] = f_odeDef_integrateAndFireModel_threeByThreeClassifier(InTron_0,dInTron_0,InTronCrit,alpha,Ispd,tau_integrate,tau_refractory,jPh0Vec,W,tVec);
        for mm = 1:3
            if length(spikeTimes{mm}(:)) > 1
                numSpikes = length(spikeTimes{mm}(2:end));%throw away the first dummy spike in the distant past
                diffSpikeTimes = diff(spikeTimes{mm}(2:end));%throw away the first dummy spike in the distant past
                averageInterSpikeInterval = sum(diffSpikeTimes)/numSpikes;
                averageSpikeRate = 1/averageInterSpikeInterval;
                outputRateBook(mm,ii,jj) = averageSpikeRate;
            else
                outputRateBook(mm,ii,jj) = 0;
            end
        end        
    end
end

givenAnswers = zeros(1,30);
for ii = 1:length(letterList)
    for jj = 1:length(defectVec)
        [~,ind] = max(outputRateBook(:,ii,jj));
        givenAnswers((ii-1)*length(defectVec)+jj) = ind;
    end
end

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(1:30,correctAnswers,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(1:30,givenAnswers,'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
ylabel('Output neuron','FontSize',fontSize,'FontName','Times')
xlabel('Input image','FontSize',fontSize,'FontName','Times')
% title(sprintf('InTron/InTronCrit for output neuron %g',jj),'FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)


return
%% analyze spike times
numSpikes = zeros(3,1);
averageSpikeRate = numSpikes;
averageInterSpikeInterval = averageSpikeRate;
for ii = 1:3
    numSpikes(ii) = length(spikeTimes{ii}(2:end));%throw away the first dummy spike in the distant past
    diffSpikeTimes = diff(spikeTimes{ii}(2:end));%throw away the first dummy spike in the distant past
    averageInterSpikeInterval(ii) = sum(diffSpikeTimes)/numSpikes(ii);
    averageSpikeRate(ii) = 1/averageInterSpikeInterval(ii);
end

fprintf(['\n\nThe letter ',letter,' was input.\nOutput neuron 1 gave %g spikes\nOutput neuron 2 gave %g spikes\nOutput neuron 3 gave %g spikes\n\n'],numSpikes(1),numSpikes(2),numSpikes(3))

%% plot InTron
for jj = 1:3
    
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    plot(tVec*1e6,InTron(jj,:)/InTronCrit,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
    ylabel('InTron/InTronCrit','FontSize',fontSize,'FontName','Times')
    xlabel('Time [\mu s]','FontSize',fontSize,'FontName','Times')
    for ii = 2:length(spikeTimes{jj})
        line([spikeTimes{jj}(ii) spikeTimes{jj}(ii)]*1e6,[-0.1 1.1],'LineStyle',':','Color','k','LineWidth',2)
    end
    ylim([-0.1 1.1])
    title(sprintf('InTron/InTronCrit for output neuron %g',jj),'FontSize',fontSize,'FontName','Times')
    set(gca,'FontSize',fontSize,'FontName',fontName)
    
end

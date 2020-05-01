%% initialize
clc
clear all
close all

[fontName,fontSize,bRGY,scrsz] = f_plotting;
fontSize = 40;

p = f_physicalConstants;

%% inputs

%parameters
Ispd = 1e-6;%current through a single snspd of the receiver array
alpha = 1;%efficiency factor scaling sent photons to received photons
tau_refractory = 50e-9;%refractory period
InTronCrit = Ispd*100;

%input signal
tau_switchVec = [10e-6 10e-3];%two values corresponding to two integration times; in this case, this sets the entire simulation time. jPh is always on

% receiver nanowire L/r time constant
tau_integrateVec = [1e-6 1e-3];

%temporal
t0 = 0;%beginning of simulation
nT = 1e4;%number of time steps

%inital conditions
InTron_0 = 0;

%% loop over jPh

% input photon current (photons/sec)
jPhCell{1} = 1e8*(1:1:10000);%short integration time
jPhCell{2} = 1e8*(1:1:10000);

firingRateCell = cell(2,1);
for pp = 1:length(tau_integrateVec)
    tau_integrate = tau_integrateVec(pp);
    tau_switch = tau_switchVec(pp);
    
    %input signal
    switchTimes = [0 tau_switch];%[0 0.2*tau_switch 1.2*tau_switch 2.2*tau_switch 3.2*tau_switch];
    tf = switchTimes(end);%total simulation time in seconds
    tVec = linspace(t0,tf,nT);
    
    for qq = 1:length(jPhCell{pp})
        fprintf('\n\nqq = %g of %g\n',qq,length(jPhCell{pp}))
        photonValues = [jPhCell{pp}(qq) jPhCell{pp}(qq)];%vector same length as switchTimes; this number of photons per second is applied at the corresponding switch time
        
        %% jPhDrive
        [jPhDrive] = f_photonDriveDef(tVec,switchTimes,photonValues);
        
        %% call time stepper
        [InTron,spikeTimes,jPhDrive] = f_odeDef_integrateAndFireModel_refractoryRise(InTron_0,InTronCrit,alpha,Ispd,tau_integrate,tau_refractory,jPhDrive,tVec);
        spikeTimes = spikeTimes(2:end);
        firingRateCell{pp}(qq) = length(spikeTimes)/(tf-t0);
        
    end
end

%% plot firing rate
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
loglog(jPhCell{1}*1e-9,firingRateCell{1}(:)*1e-6,'Color',bRGY(3,:),'LineWidth',2)
xlabel('Input photon current [10^9 photons/s]','FontSize',fontSize,'FontName',fontName)
ylabel('Neuron firing rate [MHz]','FontSize',fontSize,'FontName',fontName)
legend('\tau_{int} = 1\mu s')
set(gca,'FontSize',fontSize,'FontName',fontName)
% ylim([0 20])

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
loglog(jPhCell{2}*1e-9,firingRateCell{2}(:)*1e-6,'Color',bRGY(8,:),'LineWidth',2)
xlabel('Input photon current [10^9 photons/s]','FontSize',fontSize,'FontName',fontName)
ylabel('Neuron firing rate [MHz]','FontSize',fontSize,'FontName',fontName)
legend('\tau_{int} = 1 ms')
set(gca,'FontSize',fontSize,'FontName',fontName)
% ylim([0 20])
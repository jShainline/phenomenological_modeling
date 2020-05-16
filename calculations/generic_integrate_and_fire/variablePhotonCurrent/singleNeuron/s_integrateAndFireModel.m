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
InTronCrit = Ispd*100;

%temporal
t0 = 0;%beginning of simulation
tf = 5e-6;%total simulation time in seconds
nT = 1e4;%number of time steps
tVec = linspace(t0,tf,nT);

%input signal
tau_switch = 1e-6;
switchTimes = [0 0.2*tau_switch 1.2*tau_switch 2.2*tau_switch 3.2*tau_switch];
photonValues = [0 jPh0 0 jPh0 0];%vector same length as switchTimes; this number of photons per second is applied at the corresponding switch time

%inital conditions
InTron_0 = 0;
dInTron_0 = 0;

%% jPhDrive
[jPhDrive] = f_photonDriveDef(tVec,switchTimes,photonValues);

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(tVec*1e6,jPhDrive,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
ylabel('jPhDrive [photons/sec]','FontSize',fontSize,'FontName','Times')
xlabel('Time [\mu s]','FontSize',fontSize,'FontName','Times')
% title(sprintf('Total energy: %g fJ\nLED energy: %g fJ\nfraction = %g\ntotalPhotons = %g, assuming efficiency = %g',totEnergy*1e15,ledEnergy*1e15,ledEnergy/totEnergy,totalPhotons,eta),'FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

%% call time stepper
[dInTron_dt,InTron,spikeTimes,jPhDrive] = f_odeDef_integrateAndFireModel(InTron_0,dInTron_0,InTronCrit,alpha,Ispd,tau_integrate,tau_refractory,jPhDrive,tVec);

%% plot InTron
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(tVec*1e6,InTron/InTronCrit,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
ylabel('InTron/InTronCrit','FontSize',fontSize,'FontName','Times')
xlabel('Time [\mu s]','FontSize',fontSize,'FontName','Times')
for ii = 1:length(spikeTimes)
    line([spikeTimes(ii) spikeTimes(ii)]*1e6,[-0.1 1.1],'LineStyle',':','Color','k','LineWidth',2)
end
ylim([-0.1 1.1])
% title(sprintf('Total energy: %g fJ\nLED energy: %g fJ\nfraction = %g\ntotalPhotons = %g, assuming efficiency = %g',totEnergy*1e15,ledEnergy*1e15,ledEnergy/totEnergy,totalPhotons,eta),'FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
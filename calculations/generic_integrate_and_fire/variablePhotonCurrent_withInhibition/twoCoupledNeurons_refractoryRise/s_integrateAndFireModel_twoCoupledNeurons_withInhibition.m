%% initialize
clc
clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;

p = f_physicalConstants;

%% input parameters

%parameters
Ispd_e = 1e-6;%current through a single snspd of the excitatory array
Ispd_i = 1e-6;%current through a single snspd of the inhibitory array
numNanowires_i = 0;
Imax_i = Ispd_i*numNanowires_i;

jPh0_e = 5e7;%photons per second into excitatory port when pulse is on
jPh0_i = 0;%photons per second into inhibitory port when pulse is on
tau_integrate_e = 1e-6;%receiver integration time
tau_integrate_i = 1e-6;%receiver integration time

tau_refractory = 50e-9;%refractory period
refractoryAffectsInhibition = 'no';
I_tTronCrit = Ispd_e*30+Imax_i;
alpha = 1;%efficiency factor scaling sent photons to received photons

%input signal
tau_switch_e = 1e-6;
switchTimes1_e = [0 0.2*tau_switch_e 5.2*tau_switch_e];%[0 0.2*tau_switch 1.2*tau_switch 2.2*tau_switch 3.2*tau_switch];
photonValues1_e = [0 jPh0_e 0];%[0 jPh0 0 jPh0 0];%vector same length as switchTimes; this number of photons per second is applied at the corresponding switch time
% tau_switch_i = 1e-6;
switchTimes2_e = [0 0.2*tau_switch_e 5.2*tau_switch_e];%[0 0.2*tau_switch 1.2*tau_switch 2.2*tau_switch 3.2*tau_switch];
photonValues2_e = [0 jPh0_e 0];%[0 jPh0 0 jPh0 0];%vector same length as switchTimes; this number of photons per second is applied at the corresponding switch time

%spike events
photonsPerSpike = 30;
spikeDuration = 1000e-9;%seconds
spikeTimeDelay = 0;%seconds
typeOfCoupling = 'excitatory';%'inhibitory', 'excitatory', or 'none'

%temporal
t0 = 0;%beginning of simulation
tf = switchTimes1_e(end)+5e-6;%0.2*tau_switch_e;%total simulation time in seconds
nT = 1e4;%number of time steps
tVec = linspace(t0,tf,nT);

%inital conditions
I_tTron_0 = 0;

%% jPhDrive
[jPh1_e] = f_photonDriveDef(tVec,switchTimes1_e,photonValues1_e);
[jPh2_e] = f_photonDriveDef(tVec,switchTimes2_e,photonValues2_e);
jPhDrive_e = {jPh1_e,jPh2_e};

% [jPh_i] = f_photonDriveDef(tVec,switchTimes_i,photonValues_i);
jPhDrive_i = {zeros(nT,1),zeros(nT,1)};%{jPh_i,jPh_i};

% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(tVec*1e6,jPhDrive_e{1},'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
% ylabel('jPhDrive [photons/sec]','FontSize',fontSize,'FontName','Times')
% xlabel('Time [\mu s]','FontSize',fontSize,'FontName','Times')
% % title(sprintf('Total energy: %g fJ\nLED energy: %g fJ\nfraction = %g\ntotalPhotons = %g, assuming efficiency = %g',totEnergy*1e15,ledEnergy*1e15,ledEnergy/totEnergy,totalPhotons,eta),'FontSize',fontSize,'FontName','Times')
% set(gca,'FontSize',fontSize,'FontName',fontName)

%% call time stepper
[I_tTron,I_tTron_e,I_tTron_i,spikeTimes,jPhDrive_e,jPhCouple_e,jPhCouple_i] = f_odeDef_integrateAndFireModel_twoCoupledNeurons_withInhibition(I_tTron_0,I_tTronCrit,Imax_i,alpha,Ispd_e,Ispd_i,tau_integrate_e,tau_integrate_i,tau_refractory,jPhDrive_e,jPhDrive_i,tVec,refractoryAffectsInhibition,photonsPerSpike,typeOfCoupling,spikeTimeDelay,spikeDuration);

%% plot excitatory photon currents
figureCaptions = {sprintf('jPh_0^e = %g photons/sec',jPh0_e),...
                  sprintf(''),...
                  sprintf('jPh_0^i = %g photons/sec',jPh0_i),...
                  sprintf(''),...
                  sprintf('photonsPerSpike = %g',photonsPerSpike),...
                  sprintf(''),...
                  sprintf(['typeOfCoupling = ',typeOfCoupling]),...
                  sprintf(''),...
                  sprintf('tau_{int}^e = %g us',tau_integrate_e*1e6),...
                  sprintf(''),...
                  sprintf('tau_{int}^i = %g us',tau_integrate_i*1e6),...
                  };
              
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
[AX,H1,H2] = plotyy(tVec*1e6,jPhDrive_e{1},tVec*1e6,jPhDrive_e{2});

set(get(AX(1),'Ylabel'),'String','jPhDrive_e(1) [photons per second]','FontSize',fontSize,'FontName','Times') 
set(get(AX(2),'Ylabel'),'String','jPhDrive_e(2) [photons per second]','FontSize',fontSize,'FontName','Times') 
set(H1,'LineStyle','-','LineWidth',3,'Color',bRGY(3,:));
set(H2,'LineStyle','-.','LineWidth',2,'Color',bRGY(8,:));
set(AX,{'ycolor'},{bRGY(3,:),;bRGY(8,:)})
set(AX,'FontSize',fontSize,'FontName',fontName)

jPhDrive_e_range = max(jPhDrive_e{1})-min(jPhDrive_e{1});
set(AX(1),'ylim',[min(jPhDrive_e{1})-jPhDrive_e_range/10 max(jPhDrive_e{1})+jPhDrive_e_range/10]);

jPhDrive_e_range = max(jPhDrive_e{2})-min(jPhDrive_e{2});
set(AX(2),'ylim',[min(jPhDrive_e{2})-jPhDrive_e_range/10 max(jPhDrive_e{2})+jPhDrive_e_range/10]);

set(AX,'xlim',[tVec(1)*1e6 tVec(end)*1e6]);

k1 = gtext(figureCaptions(1:length(figureCaptions)));
set(k1,'FontSize',fontSize_legend,'FontName','Times')

%% plot coupling photon currents
% figureCaptions = {sprintf('jPh_0^e = %g photons/sec',jPh0_e),...
%                   sprintf(''),...
%                   sprintf('jPh_0^i = %g photons/sec',jPh0_i),...
%                   sprintf(''),...
%                   sprintf('photonsPerSpike = %g',photonsPerSpike),...
%                   sprintf(''),...
%                   sprintf(['typeOfCoupling = ',typeOfCoupling]),...
%                   sprintf(''),...
%                   sprintf('tau_{int}^e = %g us',tau_integrate_e*1e6),...
%                   sprintf(''),...
%                   sprintf('tau_{int}^i = %g us',tau_integrate_i*1e6),...
%                   };
%               
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% [AX,H1,H2] = plotyy(tVec*1e6,jPhCouple_e{1},tVec*1e6,jPhCouple_e{2});
% 
% set(get(AX(1),'Ylabel'),'String','jPhCouple_e(1) [photons per second]','FontSize',fontSize,'FontName','Times') 
% set(get(AX(2),'Ylabel'),'String','jPhCouple_e(2) [photons per second]','FontSize',fontSize,'FontName','Times') 
% set(H1,'LineStyle','-','LineWidth',3,'Color',bRGY(3,:));
% set(H2,'LineStyle','-.','LineWidth',2,'Color',bRGY(8,:));
% set(AX,{'ycolor'},{bRGY(3,:),;bRGY(8,:)})
% set(AX,'FontSize',fontSize,'FontName',fontName)
% 
% jPhCouple_e_range = max(jPhCouple_e{1})-min(jPhCouple_e{1});
% set(AX(1),'ylim',[min(jPhCouple_e{1})-jPhCouple_e_range/10 max(jPhCouple_e{1})+jPhCouple_e_range/10]);
% 
% jPhCouple_e_range = max(jPhCouple_e{2})-min(jPhCouple_e{2});
% set(AX(2),'ylim',[min(jPhCouple_e{2})-jPhCouple_e_range/10 max(jPhCouple_e{2})+jPhCouple_e_range/10]);
% 
% set(AX,'xlim',[tVec(1)*1e6 tVec(end)*1e6]);
% 
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')

%% plot I_tTron
             
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
[AX,H1,H2] = plotyy(tVec*1e6,I_tTron{1}/I_tTronCrit,tVec*1e6,I_tTron{2}/I_tTronCrit);
% H3 = line(tVec*1e6,jPhDrive_e{1},'Parent',AX(2),'LineStyle','-.','LineWidth',2,'Color',bRGY(13,:));
% H4 = line(tVec*1e6,jPhDrive_e{2},'Parent',AX(2),'LineStyle','-.','LineWidth',2,'Color',bRGY(18,:));

set(get(AX(1),'Ylabel'),'String','I_{tTron}(1)/I_{tTron}^c','FontSize',fontSize,'FontName','Times') 
set(get(AX(2),'Ylabel'),'String','I_{tTron}(2)/I_{tTron}^c','FontSize',fontSize,'FontName','Times') 
set(H1,'LineStyle','-','LineWidth',3,'Color',bRGY(3,:));
set(H2,'LineStyle','-.','LineWidth',2,'Color',bRGY(8,:));
set(AX,{'ycolor'},{bRGY(3,:),;bRGY(8,:)})
set(AX,'FontSize',fontSize,'FontName',fontName)
% lgd = legend([H2;H3],'j_{ph}^e(1)','j_{ph}^e(2)');
% set(lgd,'FontSize',fontSize_legend,'FontName',fontName)

xlabel('Time [\mu s]','FontSize',fontSize,'FontName','Times')
for ii = 2:length(spikeTimes{1})
    line([spikeTimes{1}(ii) spikeTimes{1}(ii)]*1e6,[-0.1 1.1],'LineStyle',':','Color',bRGY(3,:),'LineWidth',1)
end
set(AX(1),'ylim',[-0.1 1.1]);

for ii = 2:length(spikeTimes{2})
    line([spikeTimes{2}(ii) spikeTimes{2}(ii)]*1e6,[-0.1 1.1],'LineStyle',':','Color',bRGY(8,:),'LineWidth',1)
end
set(AX(2),'ylim',[-0.1 1.1]);

% title(sprintf('Total energy: %g fJ\nLED energy: %g fJ\nfraction = %g\ntotalPhotons = %g, assuming efficiency = %g',totEnergy*1e15,ledEnergy*1e15,ledEnergy/totEnergy,totalPhotons,eta),'FontSize',fontSize,'FontName','Times')
set(AX,'FontSize',fontSize,'FontName',fontName)
set(AX,'xlim',[tVec(1)*1e6 tVec(end)*1e6]);

k1 = gtext(figureCaptions(1:length(figureCaptions)));
set(k1,'FontSize',fontSize_legend,'FontName','Times')

return

%% plot I_tTron_e and I_tTron_i
              
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
[AX,H1,H2] = plotyy(tVec*1e6,I_tTron_e,tVec*1e6,I_tTron_i);

set(get(AX(1),'Ylabel'),'String','I_{tTron}^e','FontSize',fontSize,'FontName','Times') 
set(get(AX(2),'Ylabel'),'String','I_{tTron}^i','FontSize',fontSize,'FontName','Times') 
set(H1,'LineStyle','-','LineWidth',3,'Color',bRGY(3,:));
set(H2,'LineStyle','-.','LineWidth',2,'Color',bRGY(8,:));
set(AX,{'ycolor'},{bRGY(3,:),;bRGY(8,:)})
set(AX,'FontSize',fontSize,'FontName',fontName)
xlabel('Time [\mu s]','FontSize',fontSize,'FontName','Times')
set(AX,'FontSize',fontSize,'FontName',fontName)
set(AX,'xlim',[tVec(1)*1e6 tVec(end)*1e6]);
k1 = gtext(figureCaptions(1:length(figureCaptions)));
set(k1,'FontSize',fontSize_legend,'FontName','Times')

%% plot I_tTron_i and jPhDrive_i
             
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
[AX,H1,H2] = plotyy(tVec*1e6,I_tTron_i/Imax_i,tVec*1e6,jPhDrive_i);

set(get(AX(1),'Ylabel'),'String','I_{tTron}^i/I_{tTron}^{max}','FontSize',fontSize,'FontName','Times') 
set(get(AX(2),'Ylabel'),'String','j_{ph}^i [photons/sec]','FontSize',fontSize,'FontName','Times') 
set(H1,'LineStyle','-','LineWidth',3,'Color',bRGY(3,:));
set(H2,'LineStyle','-.','LineWidth',2,'Color',bRGY(8,:));
set(AX,{'ycolor'},{bRGY(3,:),;bRGY(8,:)})
set(AX,'FontSize',fontSize,'FontName',fontName)

xlabel('Time [\mu s]','FontSize',fontSize,'FontName','Times')
set(AX,'FontSize',fontSize,'FontName',fontName)
set(AX,'xlim',[tVec(1)*1e6 tVec(end)*1e6]);

k1 = gtext(figureCaptions(1:length(figureCaptions)));
set(k1,'FontSize',fontSize_legend,'FontName','Times')

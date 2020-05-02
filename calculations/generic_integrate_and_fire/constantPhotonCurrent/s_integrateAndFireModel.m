%% initialize
clc
clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;

p = f_physicalConstants;

%% model 1: varying IntC

%inputs

% receiver nanowire L/r time constant
tauSPDVec = [100e-9 1e-6];%short integration time
% tauSPDVec = [10e-6 1e-3];%long integration time

% input photon current (photons/sec)
jPhVec = 1e7*(1:1:1000);%short integration time
% jPhVec = 1e4*(1:1:100000);%long integration time


tauRef = 50e-9;%refractory period
alpha = 1;%constant of proportionality between photon current and nanowire absorption (number between 0 and 1, accounts for detector efficiency and nanowire statistics)

nNW = 1000;%number of nanowires in receiver array

Ispd = 1e-6;%single-photon detector bias current
IntCVec = nNW*Ispd*([0.01 0.05 0.1 0.5]);%critical current of nTron

% interspike interval and firing rate
% tISI = -tau*log(1-(nNW./(2*tau*jPhVec)));%interspike interval, assuming IntC = nNW*Ispd/2
tISI = zeros(length(jPhVec),length(tauSPDVec),length(IntCVec));

for ii = 1:length(tauSPDVec)
    for jj = 1:length(IntCVec)
        for kk = 1:length(jPhVec)
            if IntCVec(jj)/(tauSPDVec(ii)*alpha*jPhVec(kk)*Ispd) < 1 && IntCVec(jj)/(tauSPDVec(ii)*alpha*jPhVec(kk)*Ispd) > 0
                tISI(kk,ii,jj) = -tauSPDVec(ii)*log(1-(IntCVec(jj)/(tauSPDVec(ii)*alpha*jPhVec(kk)*Ispd)))+tauRef;%interspike interval
            else
                tISI(kk,ii,jj) = 1e12;
            end
        end
    end
end

rISI = tISI.^-1;

% plot firing rate

% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% colorVec = [1 6 11 2 7 12 3 8 13];
% legendStr = 'legend(';
% for jj = 1:length(IntCVec)
%     for ii = 1:length(tauSPDVec)
%         plot(jPhVec*1e-9,rISI(:,ii,jj)*1e-6,'Color',bRGY(colorVec( (ii-1)*length(IntCVec)+jj ),:),'LineWidth',2)
%         legendStr = [legendStr sprintf('''I_{tT}^{c} = %g n_{nw} I_{spd}, tau_{int} = %g ns''',IntCVec(jj)/nNW/Ispd,tauSPDVec(ii)*1e9) ','];
%         hold on
%     end
% end
% legendStr = [legendStr(1:end-1),');'];
% eval(legendStr);
% xlabel('Input photon current [10^9 photons/s]','FontSize',fontSize,'FontName',fontName)
% ylabel('Neuron firing rate [MHz]','FontSize',fontSize,'FontName',fontName)
% % title(sprintf('tauSPD = %g ns, tauRef = %g, nNW = %g',tauSPDVec*1e9,nNW),'FontSize',fontSize,'FontName',fontName)
% set(gca,'FontSize',fontSize,'FontName',fontName)
% ylim([0 20])


figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
colorVec = [1 6 11 16 3 8 13 18];
legendStr = 'lgd = legend(';
for jj = 1:length(IntCVec)
    for ii = 1:length(tauSPDVec)
        [ind] = find(rISI(:,ii,jj) > 1);
        semilogx(jPhVec(ind)*1e-9,rISI(ind,ii,jj)*1e-6,'Color',bRGY(colorVec( (ii-1)*length(IntCVec)+jj ),:),'LineWidth',2)
        hold on
%         legendStr = [legendStr sprintf('''I_{tT}^{c} = %g n_{nw} I_{spd}, tau_{int} = %g us''',IntCVec(jj)/nNW/Ispd,tauSPDVec(ii)*1e6) ','];%short integration time
                legendStr = [legendStr sprintf('''I_{tT}^{c} = %g n_{nw} I_{spd}, tau_{int} = %g ms''',IntCVec(jj)/nNW/Ispd,tauSPDVec(ii)*1e3) ','];%long integration time
        hold on
    end
end
legendStr = [legendStr(1:end-1),');'];
eval(legendStr);
xlabel('Input photon current [10^9 photons/s]','FontSize',fontSize,'FontName',fontName)
ylabel('Neuron firing rate [MHz]','FontSize',fontSize,'FontName',fontName)
% title(sprintf('tauSPD = %g ns, tauRef = %g, nNW = %g',tauSPDVec*1e9,nNW),'FontSize',fontSize,'FontName',fontName)
set(gca,'FontSize',fontSize,'FontName',fontName)
lgd.FontSize = fontSize_legend;
ylim([0 20])


figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
colorVec = [1 6 11 16 3 8 13 18];
legendStr = 'lgd = legend(';
for jj = 1:length(IntCVec)
    for ii = 1:length(tauSPDVec)
        [ind] = find(rISI(:,ii,jj) > 1);
        loglog(jPhVec(ind)*1e-9,rISI(ind,ii,jj)*1e-6,'Color',bRGY(colorVec( (ii-1)*length(IntCVec)+jj ),:),'LineWidth',2)
        hold on
%         legendStr = [legendStr sprintf('''I_{tT}^{c} = %g n_{nw} I_{spd}, tau_{int} = %g us''',IntCVec(jj)/nNW/Ispd,tauSPDVec(ii)*1e6) ','];%short integration time
                legendStr = [legendStr sprintf('''I_{tT}^{c} = %g n_{nw} I_{spd}, tau_{int} = %g ms''',IntCVec(jj)/nNW/Ispd,tauSPDVec(ii)*1e3) ','];%long integration time
        hold on
    end
end
legendStr = [legendStr(1:end-1),');'];
eval(legendStr);
xlabel('Input photon current [10^9 photons/s]','FontSize',fontSize,'FontName',fontName)
ylabel('Neuron firing rate [MHz]','FontSize',fontSize,'FontName',fontName)
% title(sprintf('tauSPD = %g ns, tauRef = %g, nNW = %g',tauSPDVec*1e9,nNW),'FontSize',fontSize,'FontName',fontName)
set(gca,'FontSize',fontSize,'FontName',fontName)
lgd.FontSize = fontSize_legend;
ylim([0.1 20])


%% plot number of photons per action potential
gamma = zeros(length(jPhVec),length(tauSPDVec),length(IntCVec));
for qq = 1:length(tauSPDVec)
    for pp = 1:length(IntCVec)
        for kk = 1:length(jPhVec)
            if real(rISI(kk,qq,pp)) > 1
                gamma(kk,qq,pp) = jPhVec(kk)/real(rISI(kk,qq,pp));%number of photons per spike
            end
        end
    end
end
eta = gamma/nNW;

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
colorVec = [1 6 11 16 3 8 13 18];
legendStr = 'lgd = legend(';
for jj = 1:length(IntCVec)
    for ii = 1:length(tauSPDVec)
        [ind] = find(rISI(:,ii,jj) > 1);
        loglog(real(rISI(ind,ii,jj))*1e-6,gamma(ind,ii,jj),'Color',bRGY(colorVec( (ii-1)*length(IntCVec)+jj ),:),'LineWidth',2)
        hold on
%         legendStr = [legendStr sprintf('''I_{tT}^{c} = %g n_{nw} I_{spd}, tau_{int} = %g us''',IntCVec(jj)/nNW/Ispd,tauSPDVec(ii)*1e6) ','];%short integration time
legendStr = [legendStr sprintf('''I_{tT}^{c} = %g n_{nw} I_{spd}, tau_{int} = %g ms''',IntCVec(jj)/nNW/Ispd,tauSPDVec(ii)*1e3) ','];%long integration time
        hold on
    end
end
legendStr = [legendStr(1:end-1),');'];
eval(legendStr);
xlabel('Firing rate [MHz]','FontSize',fontSize,'FontName',fontName)
ylabel('Photons per spike','FontSize',fontSize,'FontName',fontName)
% title(sprintf('tauSPD = %g ns, tauRef = %g, nNW = %g',tauSPDVec*1e9,nNW),'FontSize',fontSize,'FontName',fontName)
set(gca,'FontSize',fontSize,'FontName',fontName)
lgd.FontSize = fontSize_legend;
% ylim([0.1 20])
% xlim([0.1 20])

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
colorVec = [1 6 11 16 3 8 13 18];
legendStr = 'lgd = legend(';
for jj = 1:length(IntCVec)
    for ii = 1:length(tauSPDVec)
        [ind] = find(rISI(:,ii,jj) > 1);
        loglog(real(rISI(ind,ii,jj))*1e-6,eta(ind,ii,jj),'Color',bRGY(colorVec( (ii-1)*length(IntCVec)+jj ),:),'LineWidth',2)
        hold on
        legendStr = [legendStr sprintf('''I_{tT}^{c} = %g n_{nw} I_{spd}, tau_{int} = %g us''',IntCVec(jj)/nNW/Ispd,tauSPDVec(ii)*1e6) ','];%short integration time
% legendStr = [legendStr sprintf('''I_{tT}^{c} = %g n_{nw} I_{spd}, tau_{int} = %g ms''',IntCVec(jj)/nNW/Ispd,tauSPDVec(ii)*1e3) ','];%long integration time
        hold on
    end
end
legendStr = [legendStr(1:end-1),');'];
eval(legendStr);
xlabel('Firing rate [MHz]','FontSize',fontSize,'FontName',fontName)
ylabel('\eta [Photons/spike/nanowire]','FontSize',fontSize,'FontName',fontName)
% title(sprintf('tauSPD = %g ns, tauRef = %g, nNW = %g',tauSPDVec*1e9,nNW),'FontSize',fontSize,'FontName',fontName)
set(gca,'FontSize',fontSize,'FontName',fontName)
lgd.FontSize = fontSize_legend;
% ylim([0.1 20])
% xlim([0.1 20])
return

%% model 2: varying Ispd

%inputs

tauSPDVec = [100e-9 1e-6 10e-6];%[50e-9 100e-9 200e-9];%receiver nanowire L/r time constant
tauRef = 50e-9;%refractory period

jPhVec = 1e4*(1:1:100000);%input photon current (photons/sec)
alpha = 1;%constant of proportionality between photon current and nanowire absorption (number between 0 and 1, accounts for detector efficiency and nanowire statistics)

nNW = 1000;%number of nanowires in receiver array

IspdC = 1e-6;%critical current of individual SNSPD
IspdVec = IspdC*[0.4 0.9];%single-photon detector bias current
IntCVec = nNW*IspdC*[0.01 0.1];%critical current of nTron

% interspike interval and firing rate
% tISI = -tau*log(1-(nNW./(2*tau*jPhVec)));%interspike interval, assuming IntC = nNW*Ispd/2
tISI = zeros(length(jPhVec),length(tauSPDVec),length(IntCVec),length(IspdVec));
eta = tISI;

for qq = 1:length(tauSPDVec)
    for pp = 1:length(IntCVec)
        for jj = 1:length(IspdVec)
            for kk = 1:length(jPhVec)
                if IntCVec(pp)/(tauSPDVec(qq)*alpha*jPhVec(kk)*IspdVec(jj)) < 1 && IntCVec(pp)/(tauSPDVec(qq)*alpha*jPhVec(kk)*IspdVec(jj)) > 0
                    tISI(kk,qq,pp,jj) = -tauSPDVec(qq)*log(1-(IntCVec(pp)/(tauSPDVec(qq)*alpha*jPhVec(kk)*IspdVec(jj))))+tauRef;%interspike interval
                else
                    tISI(kk,qq,pp,jj) = 1e12;
                end
            end
        end
    end
end

rISI = tISI.^-1;

for qq = 1:length(tauSPDVec)
    for pp = 1:length(IntCVec)
        for jj = 1:length(IspdVec)
            for kk = 1:length(jPhVec)
                eta(kk,qq,pp,jj) = jPhVec(kk)/(nNW*rISI(kk,qq,pp,jj));%number of photons per connection per spike
            end
        end
    end
end
% plot firing rate

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);

colorVec = [1 2 3 4 6 7 8 9];
legendStr = 'legend(';

counter = 1;
for qq = 1:length(tauSPDVec)
    for pp = 1:length(IntCVec)
        for jj = 1:length(IspdVec)
            plot(jPhVec*1e-9,rISI(:,qq,pp,jj)*1e-6,'Color',bRGY(colorVec( counter ),:),'LineWidth',2)
            legendStr = [legendStr sprintf('''tauSPD = %g ns, IntC = %g nNW IspdC, Ispd = %g IspdC''',tauSPDVec(qq)*1e9,IntCVec(pp)/nNW/IspdC,IspdVec(jj)/IspdC) ','];
            hold on
            counter = counter+1;
            if counter > length(colorVec); counter = 1; end
        end
    end
end
legendStr = [legendStr(1:end-1),');'];
eval(legendStr);
xlabel('Input photon current [10^9 photons/s]','FontSize',fontSize,'FontName',fontName)
ylabel('Neuron firing rate [MHz]','FontSize',fontSize,'FontName',fontName)
title(sprintf('tauRef = %g ns, nNW = %g',tauRef*1e9,nNW),'FontSize',fontSize,'FontName',fontName)
set(gca,'FontSize',fontSize,'FontName',fontName)
ylim([0 20])



figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);

colorVec = [1 2 3 4 6 7 8 9];
legendStr = 'legend(';

counter = 1;
for qq = 1:length(tauSPDVec)
    for pp = 1:length(IntCVec)
        for jj = 1:length(IspdVec)
            plot(rISI(:,qq,pp,jj)*1e-6,eta(:,qq,pp,jj),'Color',bRGY(colorVec( counter ),:),'LineWidth',2)
            legendStr = [legendStr sprintf('''tauSPD = %g ns, IntC = %g nNW IspdC, Ispd = %g IspdC''',tauSPDVec(qq)*1e9,IntCVec(pp)/nNW/IspdC,IspdVec(jj)/IspdC) ','];
            hold on
            counter = counter+1;
            if counter > length(colorVec); counter = 1; end
        end
    end
end
legendStr = [legendStr(1:end-1),');'];
eval(legendStr);
xlabel('Neuron firing rate [MHz]','FontSize',fontSize,'FontName',fontName)
ylabel('Photons per connections per spike','FontSize',fontSize,'FontName',fontName)
title(sprintf('tauRef = %g ns, nNW = %g',tauRef*1e9,nNW),'FontSize',fontSize,'FontName',fontName)
set(gca,'FontSize',fontSize,'FontName',fontName)
ylim([0 20])



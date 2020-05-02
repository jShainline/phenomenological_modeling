%% initialize
clc
clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;

p = f_physicalConstants;

%%
tauRef = 1e-6;%1e-9;%refractory period, units of seconds
tauInt = 1e-6;%50e-6;%integration time, units of seconds

Ispd = 1e-6;
synapticWeight = 1;
IPulse = synapticWeight*Ispd;
NyTc = 10;
IyTc = NyTc*Ispd;

rInVec = linspace(0,100e6,1000);
% tic
rOutVec = zeros(size(rInVec));
for ii = 1:length(rInVec)
    tN1 = 1-IyTc/(tauInt*rInVec(ii)*IPulse);
    if tN1 > 0; tN2 = tN1; else tN2 = 0; end
%     tN3
    rOutVec(ii) = (tauRef-tauInt*log(tN2))^(-1);
end
% toc

% tic
% rOutVec =  ( tauRef-tauInt*real( log(1-IyTc./(tauInt*rInVec*IPulse)) ) ).^(-1);
% % rOutVec = real( ( tauRef-tauInt*log(tauInt*rInVec*IPulse-IyTc) ).^(-1) );
% rOutVec(find(1-IyTc./(tauInt*rInVec*IPulse) < 0)) = 0; 
% toc
% return

%%
figureCaptions = {sprintf('tauRef = %g us',tauRef*1e6),...
                  sprintf('tauInt = %g us',tauInt*1e6),...
                  sprintf('Ispd = %g uA',Ispd*1e6),...
                  sprintf('synapticWeight = %g',synapticWeight),...
                  sprintf('NyTc = %g',NyTc)};

              
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(rInVec*1e-6,1-IyTc./(tauInt*rInVec*IPulse),'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
% % line([times(1) times(end)]*1e9,[Tt Tt],'LineStyle','-.','LineWidth',1,'Color',bRGY(1,:))
% % line([switchTimes(2) switchTimes(2)]*1e9-switchTimes(2)*1e9,[Tg max(Temps(:,1))],'LineStyle','-.','LineWidth',1,'Color',bRGY(11,:))
% % line([switchTimes(3) switchTimes(3)]*1e9-switchTimes(2)*1e9,[Tg max(Temps(:,1))],'LineStyle','-.','LineWidth',1,'Color',bRGY(6,:))
% xlabel('Input rate [MHz]','FontSize',fontSize,'FontName','Times')
% ylabel('1-IyTc./(tauInt*rInVec*IPulse)','FontSize',fontSize,'FontName','Times')
% set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')
%               
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(rInVec*1e-6,IyTc,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
% hold on
% plot(rInVec*1e-6,tauInt*rInVec*IPulse,'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
% % line([times(1) times(end)]*1e9,[Tt Tt],'LineStyle','-.','LineWidth',1,'Color',bRGY(1,:))
% % line([switchTimes(2) switchTimes(2)]*1e9-switchTimes(2)*1e9,[Tg max(Temps(:,1))],'LineStyle','-.','LineWidth',1,'Color',bRGY(11,:))
% % line([switchTimes(3) switchTimes(3)]*1e9-switchTimes(2)*1e9,[Tg max(Temps(:,1))],'LineStyle','-.','LineWidth',1,'Color',bRGY(6,:))
% lgd = legend('IyTc','tauInt x rInVec x IPulse');
% set(lgd,'FontSize',fontSize_legend);
% xlabel('Input rate [MHz]','FontSize',fontSize,'FontName','Times')
% ylabel('1-IyTc./(tauInt*rInVec*IPulse)','FontSize',fontSize,'FontName','Times')
% set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')
% 
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(rInVec*1e-6,tauRef,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
% hold on
% plot(rInVec*1e-6,real(-tauInt*log(1-IyTc./(tauInt*rInVec*IPulse))),'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
% % line([times(1) times(end)]*1e9,[Tt Tt],'LineStyle','-.','LineWidth',1,'Color',bRGY(1,:))
% % line([switchTimes(2) switchTimes(2)]*1e9-switchTimes(2)*1e9,[Tg max(Temps(:,1))],'LineStyle','-.','LineWidth',1,'Color',bRGY(11,:))
% % line([switchTimes(3) switchTimes(3)]*1e9-switchTimes(2)*1e9,[Tg max(Temps(:,1))],'LineStyle','-.','LineWidth',1,'Color',bRGY(6,:))
% lgd = legend('tauRef','tauInt');
% set(lgd,'FontSize',fontSize_legend);
% xlabel('Input rate [MHz]','FontSize',fontSize,'FontName','Times')
% ylabel('Value','FontSize',fontSize,'FontName','Times')
% set(gca,'FontSize',fontSize,'FontName',fontName)
% % ylim([0 1])
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(rInVec*1e-6,rOutVec*1e-6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
% line([times(1) times(end)]*1e9,[Tt Tt],'LineStyle','-.','LineWidth',1,'Color',bRGY(1,:))
% line([switchTimes(2) switchTimes(2)]*1e9-switchTimes(2)*1e9,[Tg max(Temps(:,1))],'LineStyle','-.','LineWidth',1,'Color',bRGY(11,:))
% line([switchTimes(3) switchTimes(3)]*1e9-switchTimes(2)*1e9,[Tg max(Temps(:,1))],'LineStyle','-.','LineWidth',1,'Color',bRGY(6,:))
xlabel('Input rate [MHz]','FontSize',fontSize,'FontName','Times')
ylabel('Output rate [MHz]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
k1 = gtext(figureCaptions(1:length(figureCaptions)));
set(k1,'FontSize',fontSize_legend,'FontName','Times')
% xlim([times(1) times(end)]*1e9)
%% initialize
clc; close all;
% clear all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;

p = f_physicalConstants;

%% inputs

I_b2_vec = 30*1e-6;

I_sy_vec = (33.75:0.25:42)*1e-6;
         
dPhi_si_vec = [19 44 76 107 145 183 227 271 314 359 409 459 509 566 622 685 748 811 874 949 1018 1093 1169 1251 1338 1426 1521 1621 1722 1829 1942 2061 2187 2325];
n_sfq_vec = floor(dPhi_si_vec/(2*pi));    
% dI_si_cell = {[0 2.1 8.3 18.6 29.0 39.3 66.2 97.2 132.3 171.6 217.1 268.8 326.7 395.0 475.6 572.8]};


%% fit
dense_I_sy_vec = linspace(I_sy_vec(1),I_sy_vec(end),1000);
[n_sfq_fit] = polyfit(I_sy_vec,n_sfq_vec,2);
n_sfq_fit_dimensionless = polyfit((n_sfq_vec-n_sfq_vec(1))/n_sfq_vec(end),(I_sy_vec-I_sy_vec(1))/I_sy_vec(end),2);
% [n_sfq_fit_power_law,~] = polyfit(log10(n_sfq_vec),log10(I_sy_vec),1);
% gamma = n_sfq_fit_power_law(1);
% A = 10^(n_sfq_fit_power_law(2));

%%
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(I_sy_vec*1e6,n_sfq_vec,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2,'Marker','o','MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(4,:))
hold on
plot(dense_I_sy_vec*1e6,polyval(n_sfq_fit,dense_I_sy_vec),'Color',bRGY(18,:),'LineStyle','-','LineWidth',2)
% plot(dense_I_sy_vec*1e6,polyval(n_sfq_fit_dimensionless,dense_I_sy_vec),'Color',bRGY(18,:),'LineStyle','-','LineWidth',2)
% plot(dense_I_sy_vec*1e6,A*dense_I_sy_vec.^gamma,'Color',bRGY(18,:),'LineStyle','-','LineWidth',2)
xlabel('I_{sy} [\mu A]','FontSize',fontSize,'FontName','Times')
ylabel('n_{sfq} [#]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
% lgd = legend('Data',sprintf('fit, gamma = %1.2f',gamma));
lgd = legend('Data',sprintf('fit, %g x^2 + %g x + %g',n_sfq_fit(1),n_sfq_fit(2),n_sfq_fit(3)));
% lgd = legend('Data',sprintf('fit, %g x^3 + %g x^2 + %g x + %g',n_sfq_fit(1),n_sfq_fit(2),n_sfq_fit(3),n_sfq_fit(4)));
% lgd = legend('Data',sprintf('fit, %g x^2 + %g x + %g',n_sfq_fit_dimensionless(1),n_sfq_fit_dimensionless(2),n_sfq_fit_dimensionless(3)));
% lgd = legend('Data',sprintf('fit, %g x + %g',n_sfq_fit_dimensionless(1),n_sfq_fit_dimensionless(2)));
set(lgd,'FontSize',fontSize_legend,'FontName',fontName)
% lgd = legend('Ib = 25','Ib = 30','Ib = 35','fit to Ib = 30');
% set(lgd,'FontSize',fontSize_legend,'FontName',fontName)
% title_str = ['Quadratic fit in dimensionless units: ' num2str(dI_si_fit_dimensionless(1)) 'x^2 + ' num2str(dI_si_fit_dimensionless(2)) 'x + ' num2str(dI_si_fit_dimensionless(3))];
% title_str = ['Qubic fit in dimensionless units: ' num2str(n_sfq_fit_dimensionless(1)) 'x^3 + ' num2str(n_sfq_fit_dimensionless(2)) 'x^2 + ' num2str(n_sfq_fit_dimensionless(3)) 'x + ' num2str(n_sfq_fit_dimensionless(4))];
% title_str = ['Qubic fit in dimensionless units: ' num2str(n_sfq_fit_dimensionless(1)) 'x^2 + ' num2str(n_sfq_fit_dimensionless(2)) 'x + ' num2str(n_sfq_fit_dimensionless(3))];
% title(title_str,'FontSize',fontSize,'FontName',fontName)
% ylim([0 0.7])
% xlim([35.5 42.5])
grid on

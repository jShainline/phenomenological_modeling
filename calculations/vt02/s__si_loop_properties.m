%% initialize
clc
clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%% inputs
beta_L_vec = logspace(3,6,100);
I_sy_vec = linspace(27e-6,40e-6,100);
L_si_sparse_vec = [10e-9 20e-9 40e-9 80e-9 100e-9 200e-9 400e-9 800e-9 1e-6 2e-6 4e-6 8e-6 10e-6 20e-6 40e-6 80e-6];
I_b2 = 30e-6;

%% constants
Ic_jj = 40e-6;%critical current of JJ
L_si_vec = p.Phi0*beta_L_vec/(2*pi*Ic_jj);%SI loop inductance
I_si_sat = f__Isisat_vs_Ib2(I_b2);

%% plot L_si vs beta_L
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
loglog(beta_L_vec,L_si_vec*1e9,'Color',bRGY(3,:),'LineStyle','-','LineWidth',2)
ylabel('L_{si} [nH]','FontSize',fontSize,'FontName','Times')
xlabel('\beta_L','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
info_str = sprintf('I_c^{(jj)} = %g uA',Ic_jj*1e6);
title(info_str,'FontSize',20,'FontName',fontName)
plot_name = sprintf('L_si__vs__beta_L__Ic%guA.png',Ic_jj*1e6);
grid on
saveas(gcf,plot_name,'png')

%% I_si vs I_sy

I_si_cell = cell(length(L_si_sparse_vec),1);
n_fq_cell = cell(length(L_si_sparse_vec),1);
for ii = 1:length(L_si_sparse_vec)
    [n_fq_vec,I_0_a] = f__I0_vs_Isy(I_sy_vec,L_si_sparse_vec(ii));
    n_fq_cell{ii} = n_fq_vec;
    I_si_cell{ii} = min(n_fq_vec*p.Phi0/L_si_sparse_vec(ii),I_si_sat);
end

%% plot n_fq
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
lgd_str = 'lgd = legend(';
color_map = [1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19];
for pp = 1:length(L_si_sparse_vec)
    semilogy(I_sy_vec*1e6,n_fq_cell{pp},'Color',bRGY(color_map(pp),:),'LineStyle','-','LineWidth',2)
    hold on
    lgd_str = [lgd_str '''' sprintf('L_{si} = %g nH',L_si_sparse_vec(pp)*1e9) '''' ','];
end
lgd_str = [lgd_str(1:end-1) ');'];
eval(lgd_str)
lgd.FontSize = 14;
ylabel('n_{fq} [#]','FontSize',fontSize,'FontName','Times')
xlabel('I_{sy} [uA]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
info_str = sprintf('I_c^{(jj)} = %g uA',Ic_jj*1e6);
title(info_str,'FontSize',20,'FontName',fontName)
plot_name = sprintf('n_fq__vs__I_sy.png');
grid on
saveas(gcf,plot_name,'png')

%% plot Delta I_si
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
lgd_str = 'lgd = legend(';
color_map = [1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19];
for pp = 1:length(L_si_sparse_vec)
    semilogy(I_sy_vec*1e6,I_si_cell{pp}*1e6,'Color',bRGY(color_map(pp),:),'LineStyle','-','LineWidth',2)
    hold on
    lgd_str = [lgd_str '''' sprintf('L_{si} = %g nH',L_si_sparse_vec(pp)*1e9) '''' ','];
end
lgd_str = [lgd_str(1:end-1) ');'];
eval(lgd_str)
lgd.FontSize = 14;
ylabel('\Delta I_{si} [\mu A]','FontSize',fontSize,'FontName','Times')
xlabel('I_{sy} [uA]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
info_str = sprintf('I_c^{(jj)} = %g uA',Ic_jj*1e6);
title(info_str,'FontSize',20,'FontName',fontName)
plot_name = sprintf('I_si__vs__I_sy.png');
grid on
saveas(gcf,plot_name,'png')

%% norm to I_si_sat
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
lgd_str = 'lgd = legend(';
color_map = [1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19];
for pp = 1:length(L_si_sparse_vec)
    semilogy(I_sy_vec*1e6,I_si_cell{pp}/I_si_sat,'Color',bRGY(color_map(pp),:),'LineStyle','-','LineWidth',2)
    hold on
    lgd_str = [lgd_str '''' sprintf('L_{si} = %g nH',L_si_sparse_vec(pp)*1e9) '''' ','];
end
lgd_str = [lgd_str(1:end-1) ');'];
eval(lgd_str)
lgd.FontSize = 14;
ylabel('\Delta I_{si}/I_{si}^{sat}','FontSize',fontSize,'FontName','Times')
xlabel('I_{sy} [uA]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
info_str = sprintf('I_c^{(jj)} = %g uA; I_{b2} = %g uA',Ic_jj*1e6,I_b2*1e6);
title(info_str,'FontSize',20,'FontName',fontName)
plot_name = sprintf('I_si__vs__I_sy__norm.png');
grid on
saveas(gcf,plot_name,'png')
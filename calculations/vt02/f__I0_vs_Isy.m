function [n_fq,I_0] = f__I0_vs_Isy(I_sy,L_si)

I_sy_vec = (27:1:40)*1e-6;
n_fq_vec = [28 80 143 216 299 393 497 614 742 885 1042 1217 1413 1635];
[n_fq_fit] = polyfit(I_sy_vec*1e6,n_fq_vec,2);

p = f_physicalConstants;
n_fq = floor(polyval(n_fq_fit,I_sy*1e6));
I_0 = n_fq*p.Phi0/L_si;

% dense_I_sy_vec = linspace(I_sy_vec(1)*1e6,I_sy_vec(end)*1e6,1000);
% n_fq_fit_vec = polyval(n_fq_fit,dense_I_sy_vec);
% 
% fit_str = sprintf('n_{fq} = %g*I_{sy}^2+%g*I_{sy}+%g, where n_{fq} is the number of flux quanta generated in a synapse event and I_{sy} is the synaptic bias current in uA',n_fq_fit(1),n_fq_fit(2),n_fq_fit(3));
% 
% [fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(I_sy_vec*1e6,n_fq_vec,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3,'Marker','o','MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
% hold on
% plot(dense_I_sy_vec,n_fq_fit_vec,'Color',bRGY(8,:),'LineStyle','-','LineWidth',1.5)
% ylabel('n_{fq} [#]','FontSize',fontSize,'FontName','Times')
% xlabel('I_{sy} [\mu A]','FontSize',fontSize,'FontName','Times')
% set(gca,'FontSize',fontSize,'FontName',fontName)
% title(fit_str,'FontSize',14,'FontName','Times')

end
function [I_si_sat] = f__Isisat_vs_Ib2(I_b2)

I_b2_vec = (25:1:38)*1e-6;
I_si_sat_vec = [8.89 9.88 10.89 11.88 12.86 13.88 14.86 15.88 16.87 17.88 18.87 19.88 20.87 21.88]*1e-6;
[I_si_sat_fit] = polyfit(I_b2_vec*1e6,I_si_sat_vec*1e6,1);
I_si_sat = polyval(I_si_sat_fit,I_b2*1e6)*1e-6;

% fit_str = sprintf('I_{si}^{sat} = %g*I_{b2}+%g, where I_{si}^{sat} is the maximum current that can be added to the SI loop and I_{b2} is the bias to the JJ in that loop in uA',I_si_sat_fit(1),I_si_sat_fit(2));
% 
% dense_I_b2_vec = linspace(I_b2_vec(1)*1e6,I_b2_vec(end)*1e6,1000);
% I_si_sat_fit_vec = polyval(I_si_sat_fit,dense_I_b2_vec);
% 
% [fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(I_b2_vec*1e6,I_si_sat_vec*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3,'Marker','o','MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
% hold on
% plot(dense_I_b2_vec,I_si_sat_fit_vec,'Color',bRGY(8,:),'LineStyle','-','LineWidth',1.5)
% ylabel('I_{si}^{sat} [\mu A]','FontSize',fontSize,'FontName','Times')
% xlabel('I_{b2} [\mu A]','FontSize',fontSize,'FontName','Times')
% set(gca,'FontSize',fontSize,'FontName',fontName)
% title(fit_str,'FontSize',14,'FontName','Times')

end
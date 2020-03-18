%% initialize
clc
clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%% variable inputs

% pulse train / temporal
input_rate_vec = logspace(5,7,12);
num_spikes = 40;%num spikes in pulse train
t_0 = 1e-9;%time of first spike
num_pts_per_event = 20;

%si loop
tau_si_vec = logspace(-7,-5,20);%[100e-9 500e-9 1e-6 5e-6 10e-6];
beta_L = 1e4;
I_b2 = 30e-6;

%synaptic bias current
I_sy_vec = 27e-6;%[27 35 39]*1e-6;


%% constant inputs

%jj
Ic_jj = 40e-6;%critical current of JJ
r_jj = 4.125;%normal state resistance of JJ
L_jj = p.Phi0/(2*pi*Ic_jj);%inductance of JJ
L_si = p.Phi0*beta_L/(2*pi*Ic_jj);%SI loop inductance

%parameters for phenomenological model
tau_rise =  4.19e-9;

gamma1 = 0.9;
gamma2 = 0.158;
gamma3 = 3/4;

I_si_sat = f__Isisat_vs_Ib2(I_b2);

%% call analytical form

steady_state_vec = zeros(length(I_sy_vec),length(tau_si_vec),length(input_rate_vec));
for pp = 1:length(I_sy_vec)
    I_sy = I_sy_vec(pp);
    [n_fq,I_0_a] = f__I0_vs_Isy(I_sy,L_si);
    
    for qq = 1:length(tau_si_vec)
        tau_si = tau_si_vec(qq);
        r_si = L_si/tau_si;%SI loop resistance
        f_0 = 1/tau_si;
        
        for rr = 1:length(input_rate_vec)
            input_rate = input_rate_vec(rr);
            Delta_t = 1/input_rate;
            t_f = t_0+Delta_t*num_spikes;
            input_spike_times = t_0:Delta_t:t_f;
            dt = Delta_t/num_pts_per_event;
            time_vec = 0:dt:t_f;
            
            fprintf('rr = %g of %g ...\n\n',rr,length(input_rate_vec))
            
            I_si = f__synaptic_response_function(time_vec,input_spike_times,I_0_a,I_si_sat,gamma1,gamma2,gamma3,tau_rise,tau_si);
            steady_state_vec(pp,qq,rr) = sum(I_si(end-num_pts_per_event:end))/num_pts_per_event;
            
        end
    end
end

%% plot steady state
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
lgd_str = 'lgd = legend(';
color_map = 1:20;
for pp = 1:length(I_sy_vec)
    I_sy = I_sy_vec(pp);
    for qq = 1:length(tau_si_vec)
        tau_si = tau_si_vec(qq);
        ssv = zeros(length(input_rate_vec),1);
        ssv(:) = steady_state_vec(pp,qq,:);
        semilogx(input_rate_vec*1e-6,ssv*1e6,'Color',bRGY(color_map((pp-1)*length(tau_si_vec)+qq),:),'LineStyle','-','LineWidth',2,'Marker','s','MarkerFace',bRGY(color_map((pp-1)*length(tau_si_vec)+qq),:),'MarkerEdge',bRGY(color_map((pp-1)*length(tau_si_vec)+qq),:))
        hold on
        lgd_str = [lgd_str '''' sprintf('tau_{si} = %g us',tau_si*1e6) '''' ','];
    end
end
lgd_str = [lgd_str(1:end-1) ');'];
eval(lgd_str)
lgd.FontSize = 12;
ylabel('Time-averaged current [\mu A]','FontSize',fontSize,'FontName','Times')
xlabel('Input Rate [MHz]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
grid on
info_str = sprintf('I_{sy} = %g uA; beta_L = %g; beta_L/2pi = %g; L_{si} = %g nH; I_{si}^{sat} = %g uA',I_sy*1e6,beta_L,beta_L/(2*pi),L_si*1e9,I_si_sat*1e6);
title(info_str,'FontSize',20,'FontName',fontName)
plot_name = sprintf('vary_rate__Isy%guA_betaL%g_Lsi%gnH_Isisat%guA.png',I_sy*1e6,beta_L,L_si*1e9,I_si_sat*1e6);
saveas(gcf,plot_name,'png')
% close

%% initialize
clc
clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%% variable inputs

%si loop
beta_L_vec = [0.5e5 1e5 0.5e6];
tau_si_vec = linspace(0.5e-6,12e-6,8);%[100e-9 333e-9 1e-6 3.33e-6 10e-6];
I_b2 = 30e-6;

%synaptic bias current
I_sy_vec = [27 35 39]*1e-6;

% pulse train / temporal
num_spikes = 20;%num spikes in pulse train
t_0 = 1e-9;%time of first spike
input_rate = 1e6;
Delta_t = 1/input_rate;
t_f = t_0+Delta_t*num_spikes;
input_spike_times = t_0:Delta_t:t_f;
dt = Delta_t/10;
time_vec = 0:dt:t_f;

%% constant inputs

%jj
Ic_jj = 40e-6;%critical current of JJ
r_jj = 4.125;%normal state resistance of JJ
L_jj = p.Phi0/(2*pi*Ic_jj);%inductance of JJ

%parameters for phenomenological model
tau_rise =  4.19e-9;

gamma1 = 0.9;
gamma2 = 0.158;
gamma3 = 3/4;

I_si_sat = f__Isisat_vs_Ib2(I_b2);

%% call analytical form

for pp = 1:length(beta_L_vec)
    
    fprintf('\n\npp = %g of %g ...\n\n',pp,length(beta_L_vec))
    
    beta_L = beta_L_vec(pp);
    L_si = p.Phi0*beta_L/(2*pi*Ic_jj);%SI loop inductance
    
    for qq = 1:length(I_sy_vec)
        
        fprintf('\nqq = %g of %g ...\n\n',qq,length(I_sy_vec))
        
        I_sy = I_sy_vec(qq);
        [n_fq,I_0_a] = f__I0_vs_Isy(I_sy,L_si);
        
        I_si_cell = cell(length(tau_si_vec),1);
        for rr = 1:length(tau_si_vec)
            
            fprintf('rr = %g of %g ...\n\n',rr,length(tau_si_vec))
            
            tau_si = tau_si_vec(rr);
            r_si = L_si/tau_si;%SI loop resistance            
            f_0 = 1/tau_si;            
            I_si_cell{rr} = f__synaptic_response_function(time_vec,input_spike_times,I_0_a,I_si_sat,gamma1,gamma2,gamma3,tau_rise,tau_si);
            
        end
        
        %plot temporal response
        figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
        lgd_str = 'lgd = legend(';
        color_map = [2 3 7 8 12 13 17 18];
        for ii = 1:length(tau_si_vec)
            plot(time_vec*1e6,I_si_cell{ii}*1e6,'Color',bRGY(color_map(ii),:),'LineStyle','-','LineWidth',3)
            hold on
            lgd_str = [lgd_str '''' sprintf('tau_{si} = %g us',tau_si_vec(ii)*1e6) '''' ','];
        end
        lgd_str = [lgd_str(1:end-1) ');'];
        eval(lgd_str)
        lgd.FontSize = fontSize_legend;
        ylabel('Current [\mu A]','FontSize',fontSize,'FontName','Times')
        xlabel('Time [ns]','FontSize',fontSize,'FontName','Times')
        set(gca,'FontSize',fontSize,'FontName',fontName)
        info_str = sprintf('I_{sy} = %g uA; beta_L = %g; beta_L/2pi = %g; L_{si} = %g nH; n_{fq} = %g; I_{0a} = %g uA; I_{si}^{sat} = %g uA',I_sy*1e6,beta_L,beta_L/(2*pi),L_si*1e9,n_fq,I_0_a*1e6,I_si_sat*1e6);
        title(info_str,'FontSize',14,'FontName',fontName)
        plot_name = sprintf('rate__Isy_%guA_betaL%g_Lsi%gnH_nFq%g_I0a%guA_Isisat%guA.png',I_sy*1e6,beta_L,L_si*1e9,n_fq,I_0_a*1e6,I_si_sat*1e6);
        saveas(gcf,plot_name,'png')
%         close
        
    end
end


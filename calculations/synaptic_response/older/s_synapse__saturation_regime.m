%% initialize
clc
% clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%% inputs

%spd
I_0 = 10e-6;%amplitude of current in spd pulse
tau_spd = 50e-9;

%jj
Ic_jj = 40e-6;%critical current of JJ
r_jj = 6.25;%normal state resistance of JJ
L_jj = 8.2e-12;%inductance of JJ

%si loop
I_b_0 = 35e-6;%dc bias to junction
L_si = 100e-9;%SI loop inductance
r_si = 2;%SI loop resistance

input_spike_times = [10 60 110 160 210 260]*1e-9;%times at which input synptic spikes occur

total_sim_time = input_spike_times(end)+2*tau_spd;
num_time_steps = 1000;
time_vec = linspace(0,total_sim_time,num_time_steps);%[0 totalSimTime];
initialConditions = 0;%initial current in the SI loop;

% [dI_si_dt] = f_odeDef_synapticLeakyIntegrator(t,I_si,L_si,Ic_jj,r_jj,L_jj,input_spike_times,I_0,tau_spd,p);

%% call drive function
[I_b] = f_synapse_drive_def(time_vec,input_spike_times,I_b_0,I_0,tau_spd);  
% 
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% plot(t_span*1e9,I_b*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)

%% call ode
fprintf('calling odeX...\n\n')
pause(0.2)
tic
[time_vec,I_si_vec] = ode45(@(t,I_si_vec) f__ode_def__synaptic_leaky_integrator(t,I_si_vec,L_si,r_si,Ic_jj,r_jj,L_jj,input_spike_times,I_b_0,I_0,tau_spd,p),time_vec,initialConditions,odeset('RelTol',1e-12,'AbsTol',1e-12));
toc
fprintf('\n\nit took %g seconds to run ode45\n\n',toc)

%% call analytical form
tau_rise = 0.8e-9;
tau_fall = L_si/r_si;%(L_jj+L_si)/r_si;
I_0_a = 3.5*max(max(I_si_vec));
I_si_vec_analytical = f__synaptic_response_function__saturation_regime(time_vec,input_spike_times,I_0_a,tau_rise,tau_fall);

%% read wrSpice
read_wr_data = 'yes';
if strcmp(read_wr_data,'yes')
    
    fprintf('\nreading WRSpice data...\n')
    file_path = 'wrSpice_data/';
    file_name = '_dat__sffg_3jj__saturation_regime_pulse_train__si10nH_50ns';
    fileID = fopen([file_path file_name],'r');
    C = textscan(fileID,'%s');
    A = C{1};
    numRows = length(A);
    tN = find(strcmp(A,'Variables:'));
    tN2 = tN(2);
    numVars = str2double(A{tN(1)+1});
    
    tN = find(strcmp(A,'Points:'));
    numPts = str2double(A{tN(1)+1});
    
    tN = tN2;
    varList = cell(numVars);
    for ii = 1:numVars
        varList{ii} = A{tN2+2+3*(ii-1)};
    end
    
    tN = find(strcmp(A,'Values:'));
    data_mat = zeros(numPts,numVars);
    tN = tN+2;
    for ii = 1:numPts
        for jj = 1:numVars
            data_mat(ii,jj) = str2double(A{tN+jj-1});
        end
        tN = tN+1+numVars;
    end
    ST = fclose(fileID);
    fprintf('\ndone reading WRSpice data...\n')
    
end

%% plot
figureCaptions = {sprintf('SI loop'),...
                  };
              
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(time_vec*1e9,I_si_vec_analytical*1e6,'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)%,'Marker','o','MarkerFaceColor',bRGY(5,:),'MarkerEdgeColor','k','MarkerSize',2)
hold on
% plot(time_vec*1e9,I_si_vec*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)%,'Marker','o','MarkerFaceColor',bRGY(5,:),'MarkerEdgeColor','k','MarkerSize',2)
% plot(time_vec*1e9,(I_b-I_b_0)*1e6,'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)%,'Marker','o','MarkerFaceColor',bRGY(5,:),'MarkerEdgeColor','k','MarkerSize',2)
plot(data_mat(:,1)*1e9,data_mat(:,6)*1e6,'Color',bRGY(18,:),'LineStyle','-','LineWidth',3)
for ii = 1:length(input_spike_times)
    line([input_spike_times(ii) input_spike_times(ii)]*1e9,[min(I_si_vec) max(I_si_vec)]*1e6,'LineStyle','-.','LineWidth',1,'Color',bRGY(1,:))
end
% line([time_vec(1) time_vec(end)]*1e9,[I_0_a I_0_a]*1e6)
% line([times(1) times(end)]*1e9,[Tt Tt],'LineStyle','-.','LineWidth',1,'Color',bRGY(1,:))
% lgd = legend('I_{si}','I_b','I_{si_a}','WR');
lgd = legend('I_{si_a}','wr');
lgd.FontSize = fontSize_legend;
ylabel('Current [\mu A]','FontSize',fontSize,'FontName','Times')
xlabel('Time [ns]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')
% xlim([times(1) times(end)]*1e9)
% xlim([-1 15])
% ylim([4 7])

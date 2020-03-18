%% initialize
clc
% clear all
% close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%% inputs

%spd
I_0 = 10e-6;%amplitude of current in spd pulse
tau_spd = 10e-9;

%jj
Ic_jj = 40e-6;%critical current of JJ
r_jj = 6.25;%normal state resistance of JJ
L_jj = 8.2e-12;%inductance of JJ

%si loop
I_b_0 = 36e-6;%dc bias to junction
L_si = 50e-9;%SI loop inductance
r_si = 0.1;%SI loop resistance

input_spike_times = [10:25:302 510:25:790]*1e-9;%times at which input synptic spikes occur

total_sim_time = 1e-6;%input_spike_times(end)+5*tau_spd;
num_time_steps = 10000;
time_vec = linspace(0,total_sim_time,num_time_steps);%[0 totalSimTime];
initialConditions = 0;%initial current in the SI loop;


%% call drive function
[I_b] = f_synapse_drive_def(time_vec,input_spike_times,I_b_0,I_0,tau_spd);  
% 
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(time_vec*1e9,I_b*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)

%% call ode
fprintf('calling odeX...\n\n')
pause(0.2)
tic
[time_vec,I_si_vec] = ode45(@(t,I_si_vec) f__ode_def__synaptic_leaky_integrator(t,I_si_vec,L_si,r_si,Ic_jj,r_jj,L_jj,input_spike_times,I_b_0,I_0,tau_spd),time_vec,initialConditions,odeset('RelTol',1e-12,'AbsTol',1e-12));
toc
fprintf('\n\nit took %g seconds to run ode45\n\n',toc)

%% read wrSpice
read_wr_data = 'no';
if strcmp(read_wr_data,'yes')
    
    fprintf('\nreading WRSpice data...\n')
    file_path = 'wrSpice_data/';
%     file_name = '_dat__sffg_3jj__nonlinear_regime_pulse_train__si25nH_500ns';
file_name = '_dat__sffg_3jj__nonlinear_regime_pulse_train__si50nH_500ns';
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

%% call analytical form
I_si_wr = data_mat(:,6);
% [peaks,ind_peaks] = findpeaks(I_si_wr);
tau_rise = 1.1e-9;
tau_fall = L_si/r_si;%(L_jj+L_si)/r_si;
I_0_a = 0.21*max(max(I_si_wr));%peaks(1);%
I_si_sat = max(I_si_wr);
gamma1 = 2;%f(I_si) exponent
gamma2 = 0.5;
I_si_vec_analytical = f__synaptic_response_function__nonlinear_regime(time_vec,input_spike_times,I_0_a,I_si_sat,gamma1,gamma2,tau_rise,tau_fall);

%% plot
figureCaptions = {sprintf('SI loop'),...
                  };
              
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(data_mat(:,1)*1e9,I_si_wr*1e6,'Color',bRGY(18,:),'LineStyle','-','LineWidth',3)
hold on
plot(time_vec*1e9,I_si_vec*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)%,'Marker','o','MarkerFaceColor',bRGY(5,:),'MarkerEdgeColor','k','MarkerSize',2)
% plot(time_vec*1e9,(I_b-I_b_0)*1e6,'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)%,'Marker','o','MarkerFaceColor',bRGY(5,:),'MarkerEdgeColor','k','MarkerSize',2)
plot(time_vec*1e9,I_si_vec_analytical*1e6,'Color',bRGY(13,:),'LineStyle','-','LineWidth',2)%,'Marker','o','MarkerFaceColor',bRGY(5,:),'MarkerEdgeColor','k','MarkerSize',2)
for ii = 1:length(input_spike_times)
    line([input_spike_times(ii) input_spike_times(ii)]*1e9,1.1*[min(I_si_wr) max(I_si_wr)]*1e6,'LineStyle','-.','LineWidth',1,'Color',bRGY(21,:))
end
% line([time_vec(1) time_vec(end)]*1e9,[I_0_a I_0_a]*1e6)
% line([times(1) times(end)]*1e9,[Tt Tt],'LineStyle','-.','LineWidth',1,'Color',bRGY(1,:))
% lgd = legend('I_{si}','I_b','I_{si_a}','WR');
% lgd = legend('I_{si_a}','I_{si}','wr');
lgd = legend('wr','ode','I_{si_a}');
lgd.FontSize = fontSize_legend;
ylabel('Current [\mu A]','FontSize',fontSize,'FontName','Times')
xlabel('Time [ns]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')
% xlim([times(1) times(end)]*1e9)
% xlim([-1 15])
% ylim([0 1.1*max(I_si_wr)*1e6])

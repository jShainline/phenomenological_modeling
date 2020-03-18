%% initialize
clc
% clear all
close all

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
L_si = 50e-9;%SI loop inductance
r_si = 0.5;%SI loop resistance

I_sy = 36.5e-6;

input_rate_vec = [25 50 75 100]*1e-9;

%% read wrSpice
read_wr_data = 'yes';
file_name = cell(length(input_rate_vec),1);
if strcmp(read_wr_data,'yes')
    
    fprintf('\nreading WRSpice data...\n\n')
    file_path = '../wrSpice_data/';
    
    file_name = {'_dat_sffg_3jj__Isy36p50uA_Lsi50nH_rsi0p5Ohm__25ns','_dat_sffg_3jj__Isy36p50uA_Lsi50nH_rsi0p5Ohm__50ns','_dat_sffg_3jj__Isy36p50uA_Lsi50nH_rsi0p5Ohm__75ns','_dat_sffg_3jj__Isy36p50uA_Lsi50nH_rsi0p5Ohm__100ns'};
    
    data_mat = cell(length(file_name),1);
    for qq = 1:length(file_name)
        
        fprintf('reading file %g of %g\n',qq,length(file_name))
        
        fileID = fopen([file_path file_name{qq}],'r');
        C = textscan(fileID,'%s');
        A = C{1};
        numRows = length(A);
        tN = find(strcmp(A,'Variables:'));
        tN2 = tN(2);
        numVars = str2double(A{tN(1)+1});
        
        tN = find(strcmp(A,'Points:'));
        numPts = str2double(A{tN(1)+1});
        varList = cell(numVars);
        for ii = 1:numVars
            varList{ii} = A{tN2+2+3*(ii-1)};
        end
        
        tN = find(strcmp(A,'Values:'));
        data_mat{qq} = zeros(numPts,numVars);
        tN = tN+2;
        for ii = 1:numPts
            for jj = 1:numVars
                data_mat{qq}(ii,jj) = str2double(A{tN+jj-1});
            end
            tN = tN+1+numVars;
        end
        ST = fclose(fileID);
        
    end
    fprintf('\ndone reading WRSpice data...\n')
end

%% call analytical form

input_spike_times = {(10:25:1000)*1e-9,(10:50:1000)*1e-9,(10:75:1000)*1e-9,(10:100:1000)*1e-9};%times at which input synptic spikes occur

tau_rise =  4.19e-9;
tau_fall = L_si/r_si;%(L_jj+L_si)/r_si;
I_0_a = 2.767e-6;
I_si_sat = 10e-6;%[12.2e-6 11.4e-6 11.4e-6];
gamma1_range = [0.9 0.9];%f(I_si) exponent
gamma2_range = [0.158 0.158];%f(I_si) exponent
num_gamma1 = 1;
num_gamma2 = 1;
gamma1_vec = linspace(gamma1_range(1),gamma1_range(2),num_gamma1);
gamma2_vec = linspace(gamma2_range(1),gamma2_range(2),num_gamma2);
gamma3 = 3/4;
downsampling_factor = 100;%integer determining how many time_vec points are used in fit

I_si_vec_analytical = cell(length(data_mat),1);
min_indices = cell(length(data_mat),1);
gamma1_vec_optimized = zeros(length(data_mat),1);
gamma2_vec_optimized = zeros(length(data_mat),1);

error_cell = cell(length(data_mat),1);
which_one = 1:length(data_mat);
for ii = which_one%1:length(data_mat)
    fprintf('\n\ncalculating analytical function %g of %g ...',ii,length(data_mat))
    
    time_vec = data_mat{ii}((1:downsampling_factor:end),1);
    target_vec = data_mat{ii}((1:downsampling_factor:end),2);
    error_min = 1e121;
    error_cell{ii} = zeros(num_gamma1,num_gamma2);
    for pp = 1:length(gamma1_vec)
        fprintf('\n\n iterating gamma 1, pp = %g of %g',pp,length(gamma1_vec))
        
        for qq = 1:length(gamma2_vec)
            fprintf('\n iterating gamma 2, qq = %g of %g',qq,length(gamma2_vec))
            
            I_si_vec_analytical_temp = f__synaptic_response_function__nonlinear_regime(time_vec,input_spike_times{ii},I_0_a,I_si_sat,gamma1_vec(pp),gamma2_vec(qq),gamma3,tau_rise,tau_fall);
            error = sum(abs(target_vec-I_si_vec_analytical_temp).^2);
            error_cell{ii}(pp,qq) = error;
            if error < error_min
                error_min = error;
                min_indices{ii} = [pp qq];
            end
            
        end
    end
    
    I_si_vec_analytical{ii} = f__synaptic_response_function__nonlinear_regime(time_vec,input_spike_times{ii},I_0_a,I_si_sat,gamma1_vec(min_indices{ii}(1)),gamma2_vec(min_indices{ii}(2)),gamma3,tau_rise,tau_fall);
    gamma1_vec_optimized(ii) = gamma1_vec(min_indices{ii}(1));
    gamma2_vec_optimized(ii) = gamma2_vec(min_indices{ii}(2));
    
end
fprintf('\n\n')

%% plot fits

for ii = which_one%1:length(data_mat)
    
    %     sprintf(file_name{ii}),...
    figureCaptions = {sprintf('input_rate = %2.2f uA',input_rate_vec(ii)*1e6),...
        sprintf('gamma1 = %g',gamma1_vec(min_indices{ii}(1))),...
        sprintf('gamma2 = %g',gamma2_vec(min_indices{ii}(2))),...
        };
    
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    plot(data_mat{ii}(:,1)*1e9,data_mat{ii}(:,2)*1e6,'Color',bRGY(18,:),'LineStyle','-','LineWidth',3)
    hold on
    plot(time_vec*1e9,I_si_vec_analytical{ii}*1e6,'Color',bRGY(13,:),'LineStyle','-','LineWidth',2)%,'Marker','o','MarkerFaceColor',bRGY(5,:),'MarkerEdgeColor','k','MarkerSize',2)
    for jj = 1:length(input_spike_times{ii})
        line([input_spike_times{ii}(jj) input_spike_times{ii}(jj)]*1e9,1.1*[min(data_mat{ii}(:,2)) max(data_mat{ii}(:,2))]*1e6,'LineStyle','-.','LineWidth',1,'Color',bRGY(21,:))
    end
    lgd = legend('wr','I_{si_a}');
    lgd.FontSize = fontSize_legend;
    ylabel('Current [\mu A]','FontSize',fontSize,'FontName','Times')
    xlabel('Time [ns]','FontSize',fontSize,'FontName','Times')
    set(gca,'FontSize',fontSize,'FontName',fontName)
    k1 = gtext(figureCaptions(1:length(figureCaptions)));
    set(k1,'FontSize',fontSize_legend,'FontName','Times')
    % xlim([times(1) times(end)]*1e9)
    % xlim([-1 15])
    % ylim([0 1.1*max(I_si_wr)*1e6])
    
end

%% plot errors


for ii = which_one%1:length(data_mat)
    
    %     sprintf(file_name{ii}),...
    figureCaptions = {sprintf('input_rate = %2.2f uA',input_rate_vec(ii)*1e6),...
        sprintf('gamma1 = %g',gamma1_vec(min_indices{ii}(1))),...
        sprintf('gamma2 = %g',gamma2_vec(min_indices{ii}(2))),...
        };
    
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    imagesc(gamma1_vec,gamma2_vec,log10(error_cell{ii}).')
    colorbar
    set(gca,'YDir','normal')
    xlabel('\gamma_1','FontSize',fontSize,'FontName','Times')
    ylabel('\gamma_2','FontSize',fontSize,'FontName','Times')
    title('log_{10}(error)','FontSize',fontSize,'FontName','Times')
    set(gca,'FontSize',fontSize,'FontName',fontName)
    k1 = gtext(figureCaptions(1:length(figureCaptions)));
    set(k1,'FontSize',fontSize_legend,'FontName','Times')
    % xlim([times(1) times(end)]*1e9)
    % xlim([-1 15])
    % ylim([0 1.1*max(I_si_wr)*1e6])
    
end




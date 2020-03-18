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
r_si = 1;%SI loop resistance

I_sy_vec = (34:0.5:39)*1e-6;

tau_fall = L_si/r_si;%(L_jj+L_si)/r_si;

%% read wrSpice
read_wr_data = 'no';
file_name = cell(length(I_sy_vec),1);
if strcmp(read_wr_data,'yes')
    
    fprintf('\nreading WRSpice data...\n\n')
    file_path = '../wrSpice_data/';
    
    file_name_base1 = '_dat_sffg_3jj__Isy';
    file_name_base2 = 'uA_Lsi50nH_rsi1Ohm__25ns';
    for ii = 1:length(I_sy_vec)
        temp_num = I_sy_vec(ii)*1e6;
        temp_str = num2str(temp_num,'%2.2f');
        file_name{ii} = [file_name_base1 temp_str(1:2) 'p' temp_str(4:5) file_name_base2];
    end
    
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

input_spike_times = [10:25:302 510:25:790]*1e-9;%times at which input synptic spikes occur

I_0_guess_range_cell = {[275e-9 325e-9],[0.65e-6 0.75e-6],[1.0e-6 1.2e-6],[1.5e-6 1.7e-6],[2.1e-6 2.3e-6],[2.7e-6 2.8e-6],[3.3e-6 3.4e-6],[3.9e-6 4.1e-6],[4.6e-6 4.7e-6],[5.3e-6 5.4e-6],[6e-6 6.1e-6]};
num_I_0_guess = 10;

I_si_sat_guess_range_cell = {[11.1e-6 11.1e-6],[11.1e-6 11.1e-6],[11.1e-6 11.1e-6],[11.1e-6 11.1e-6],[11.1e-6 11.1e-6],[11.1e-6 11.1e-6],[11.1e-6 11.1e-6],[11.1e-6 11.1e-6],[11.1e-6 11.1e-6],[11.1e-6 11.1e-6],[11.1e-6 11.1e-6]};
num_I_si_sat_guess = 1;

tau_rise_guess_range_cell = {[0.5e-9 1.5e-9],[1.5e-9 2e-9],[2e-9 3e-9],[2.75e-9 3.25e-9],[3e-9 4e-9],[4e-9 5e-9],[4.75e-9 5.25e-9],[5.25e-9 5.75e-9],[6e-9 6.5e-9],[7e-9 8e-9],[7e-9 8e-9]};
num_tau_rise_guess = 10;

gamma1 = 0.9;
gamma2 = 0.158;
gamma3 = 3/4;
downsampling_factor = 100;%integer determining how many time_vec points are used in fit

I_si_vec_analytical = cell(length(data_mat),1);
min_indices = cell(length(data_mat),1);
tau_rise_vec_optimized = zeros(length(data_mat),1);
I_0_vec_optimized = zeros(length(data_mat),1);
I_si_sat_vec_optimized = zeros(length(data_mat),1);

error_cell = cell(length(data_mat),1);
which_one = 1:length(data_mat);
for ii = which_one%1:length(data_mat)
    fprintf('\n\ncalculating analytical function %g of %g ...',ii,length(data_mat))
    
    time_vec = data_mat{ii}((1:downsampling_factor:end),1);
    target_vec = data_mat{ii}((1:downsampling_factor:end),2);
    
    I_0_vec = linspace(I_0_guess_range_cell{ii}(1),I_0_guess_range_cell{ii}(2),num_I_0_guess);
    I_si_sat_vec = linspace(I_si_sat_guess_range_cell{ii}(1),I_si_sat_guess_range_cell{ii}(2),num_I_si_sat_guess);
    tau_rise_vec = linspace(tau_rise_guess_range_cell{ii}(1),tau_rise_guess_range_cell{ii}(2),num_tau_rise_guess);
    
    error_min = 1e121;
    error_cell{ii} = zeros(num_I_0_guess,num_I_si_sat_guess,num_tau_rise_guess);
    for pp = 1:num_I_0_guess
        fprintf('\n\n iterating I_0, pp = %g of %g',pp,num_I_0_guess)
        
        for qq = 1:num_I_si_sat_guess
            fprintf('\n\n iterating tau_rise, qq = %g of %g\n',qq,num_I_si_sat_guess)
            
            for rr = 1:num_tau_rise_guess
                fprintf('\n iterating tau_rise, rr = %g of %g',rr,num_tau_rise_guess)
                
                I_si_vec_analytical_temp = f__synaptic_response_function__nonlinear_regime(time_vec,input_spike_times,I_0_vec(pp),I_si_sat_vec(qq),gamma1,gamma2,gamma3,tau_rise_vec(rr),tau_fall);
                dt = diff(time_vec);
                error = sum( (abs( target_vec(2:end)-I_si_vec_analytical_temp(2:end) ).^2).*dt )/sum( (abs( target_vec(2:end) ).^2).*dt );
                error_cell{ii}(pp,qq,rr) = error;
                if error < error_min
                    error_min = error;
                    min_indices{ii} = [pp qq rr];
                end
                
            end
        end
    end
    
    I_si_vec_analytical{ii} = f__synaptic_response_function__nonlinear_regime(time_vec,input_spike_times,I_0_vec(min_indices{ii}(1)),I_si_sat_vec(min_indices{ii}(2)),gamma1,gamma2,gamma3,tau_rise_vec(min_indices{ii}(3)),tau_fall);    
    
    I_0_vec_optimized(ii) = I_0_vec(min_indices{ii}(1));
    I_si_sat_vec_optimized(ii) = I_si_sat_vec(min_indices{ii}(2));
    tau_rise_vec_optimized(ii) = tau_rise_vec(min_indices{ii}(3));
    
end
fprintf('\n\n')

%% plot fits

for ii = which_one%1:length(data_mat)
    
    %     sprintf(file_name{ii}),...
    figureCaptions = {sprintf('I_{sy} = %2.2f uA',I_sy_vec(ii)*1e6),...
        sprintf('I_0 = %g uA',1e6*I_0_vec_optimized(ii)),...
        sprintf('I_{si}^{sat} = %g uA',1e6*I_si_sat_vec_optimized(ii)),...
        sprintf('tau_{rise} = %g ns',1e9*tau_rise_vec_optimized(ii)),...
        };
    
    figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
    plot(data_mat{ii}(:,1)*1e9,data_mat{ii}(:,2)*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
    hold on
    plot(time_vec*1e9,I_si_vec_analytical{ii}*1e6,'Color',bRGY(8,:),'LineStyle','-','LineWidth',1.5)%,'Marker','o','MarkerFaceColor',bRGY(5,:),'MarkerEdgeColor','k','MarkerSize',2)
    for jj = 1:length(input_spike_times)
        line([input_spike_times(jj) input_spike_times(jj)]*1e9,1.1*[min(data_mat{ii}(:,2)) max(data_mat{ii}(:,2))]*1e6,'LineStyle','-.','LineWidth',1,'Color',bRGY(21,:))
    end
    ylim(1.1*[min(data_mat{ii}(:,2)) max(data_mat{ii}(:,2))]*1e6)
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


% for ii = which_one%1:length(data_mat)
%     
%     %     sprintf(file_name{ii}),...
% %     figureCaptions = {sprintf('I_{sy} = %2.2f uA',I_sy_vec(ii)*1e6),...
% %         sprintf('gamma1 = %g',gamma1_vec(min_indices{ii}(1))),...
% %         sprintf('gamma2 = %g',gamma2_vec(min_indices{ii}(2))),...
% %         };
% 
%     figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
%     imagesc(I_0_vec_optimized,I_si_sat_vec_optimized,log10(error_cell{ii}(:,:,1)).')
%     colorbar
%     set(gca,'YDir','normal')
%     xlabel('I_0','FontSize',fontSize,'FontName','Times')
%     ylabel('I_{si}^{sat}','FontSize',fontSize,'FontName','Times')
%     title('log_{10}(error)','FontSize',fontSize,'FontName','Times')
%     set(gca,'FontSize',fontSize,'FontName',fontName)
% %     k1 = gtext(figureCaptions(1:length(figureCaptions)));
% %     set(k1,'FontSize',fontSize_legend,'FontName','Times')
%     % xlim([times(1) times(end)]*1e9)
%     % xlim([-1 15])
%     % ylim([0 1.1*max(I_si_wr)*1e6])
%     
% end

%% plot I_0 vs I_sy

dense_I_sy_vec = linspace(I_sy_vec(1)*1e6,I_sy_vec(end)*1e6,1000);
[I_0_fit] = polyfit(I_sy_vec*1e6,1e6*I_0_vec_optimized.',2);
% n_sfq_fit_dimensionless = polyfit((n_sfq_vec-n_sfq_vec(1))/n_sfq_vec(end),(I_sy_vec-I_sy_vec(1))/I_sy_vec(end),2);

I_0_fit_vec = polyval(I_0_fit,dense_I_sy_vec);
I_sy_test_case = 35.467;
[~,ind] = min( abs(dense_I_sy_vec - I_sy_test_case) );
I_0_test_case = I_0_fit_vec(ind);

I_0_test2 = polyval(I_0_fit,I_sy_test_case);

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(I_sy_vec*1e6,I_0_vec_optimized*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3,'Marker','o','MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
hold on
plot(dense_I_sy_vec,I_0_fit_vec,'Color',bRGY(8,:),'LineStyle','-','LineWidth',1.5)
lgd = legend('Simulation',sprintf('Fit: I_{0}(I_{sy}) =  %g I_{sy}^2 + %g I_{sy} + %g',I_0_fit(1),I_0_fit(2),I_0_fit(3)));
% lgd = legend('wr','I_{si_a}');
% lgd.FontSize = fontSize_legend;
ylabel('I_0 [\mu A]','FontSize',fontSize,'FontName','Times')
xlabel('I_{sy} [\mu A]','FontSize',fontSize,'FontName','Times')
% title(sprintf('At Isy = %g uA, I0 = %g',I_sy_test_case,I_0_test_case),'FontSize',fontSize,'FontName','Times')
% line([I_sy_test_case I_sy_test_case],[min(I_0_vec_optimized) max(I_0_vec_optimized)]*1e6)
% line([min(I_sy_vec) max(I_sy_vec)]*1e6,[I_0_test_case I_0_test_case])
set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')

%% plot min_error vs I_sy

error_vec = zeros(length(which_one),1);
for ii = which_one    
    error_vec(ii) = error_cell{ii}(min_indices{ii}(1),min_indices{ii}(2),min_indices{ii}(3));
end
    
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(I_sy_vec*1e6,error_vec,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3,'Marker','o','MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
% lgd = legend('wr','I_{si_a}');
% lgd.FontSize = fontSize_legend;
ylabel('Error','FontSize',fontSize,'FontName','Times')
xlabel('I_{sy} [\mu A]','FontSize',fontSize,'FontName','Times')
% title(sprintf('At Isy = %g uA, I0 = %g',I_sy_test_case,I_0_test_case),'FontSize',fontSize,'FontName','Times')
% line([I_sy_test_case I_sy_test_case],[min(I_0_vec_optimized) max(I_0_vec_optimized)]*1e6)
% line([min(I_sy_vec) max(I_sy_vec)]*1e6,[I_0_test_case I_0_test_case])
set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')

%% plot I_si_sat vs I_sy
figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(I_sy_vec*1e6,I_si_sat_vec_optimized*1e6,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3,'Marker','o','MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
% lgd = legend('wr','I_{si_a}');
% lgd.FontSize = fontSize_legend;
ylabel('I_0 [\mu A]','FontSize',fontSize,'FontName','Times')
xlabel('I_{si}^{sat} [\mu A]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')

%% plot tau_rise vs I_sy

[tau_rise_fit] = polyfit(I_sy_vec*1e6,1e9*tau_rise_vec_optimized.',1);
tau_rise_fit_vec = polyval(tau_rise_fit,dense_I_sy_vec);

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(I_sy_vec*1e6,tau_rise_vec_optimized*1e9,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3,'Marker','o','MarkerFaceColor',bRGY(1,:),'MarkerEdgeColor',bRGY(5,:))
hold on
plot(dense_I_sy_vec,tau_rise_fit_vec,'Color',bRGY(8,:),'LineStyle','-','LineWidth',1.5)
lgd = legend('Optimal',sprintf('Fit: tau_{rise}(I_{sy}) = %g I_{sy} + %g',tau_rise_fit(1),tau_rise_fit(2)));
% lgd = legend('wr','I_{si_a}');
% lgd.FontSize = fontSize_legend;
ylabel('\tau_{rise} [ns]','FontSize',fontSize,'FontName','Times')
xlabel('I_{sy} [\mu A]','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')
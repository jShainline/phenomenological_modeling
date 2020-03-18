%% initialize
clc
% clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%% read WR

read_wr_data = 'no';
if strcmp(read_wr_data,'yes')
    
    fprintf('\nreading WRSpice data...\n')
    file_path = 'wrSpice_data/';
    file_name = '_dat__sffg_1jj__direct_current_drive__investigating_rate_fq__si1uH_tauInf';
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

%% find num fluxons for each bias condition in WRSpice data
dc_drive_vec = [41 43 45]*1e-6;

times_switch = [10 15 20 25]*1e-9;

indices_switch = zeros(size(times_switch));
for ii = 1:length(indices_switch)
    [~,indices_switch(ii)] = min( abs(times_switch(ii)-data_mat(:,1)) );
end

num_fluxons_vec = zeros(size(dc_drive_vec));
for ii = 1:length(num_fluxons_vec)
    num_fluxons_vec(ii) = floor((data_mat(indices_switch(ii+1),5)-data_mat(indices_switch(ii),5))/(2*pi));
end

%% calculate num fluxons based on current-biased rate model
Ic_jj = 40e-6;
r_jj = 6.25;
r_fq = (r_jj/p.Phi0).*sqrt(dc_drive_vec.^2-Ic_jj^2);
num_fluxons_vec_rate_model = floor(r_fq*5e-9);

%% difference

fprintf('\n\nn_fq_wr - n_fq_rate = %g\n\n',num_fluxons_vec-num_fluxons_vec_rate_model)
fprintf('\n\n(n_fq_wr - n_fq_rate)/n_fq_wr = %g\n\n',(num_fluxons_vec-num_fluxons_vec_rate_model)./num_fluxons_vec)


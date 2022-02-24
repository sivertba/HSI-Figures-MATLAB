%% This script makes sets of test paramters for compress_scores(...)
clc; close all; clear;

% Script flags
print_check = false;
do_fft = false;
do_wavelet = true;
do_j2k = false;

parameters_cell = {};

%% FFT parameters
% Assuming fft2

method = 'fft';
tresh = [0.000005 0.00001 0.00005 0.0001];

for cr_ii = 1:length(tresh)
    if not(do_fft)
        break;
    end
    
    % Set parameter struct
    parameter_struct = struct;
    parameter_struct.method = method;
    parameter_struct.tresh = tresh(cr_ii);
    
    % Store parameter struct
    parameters_cell{end+1} = parameter_struct;
end


%% Wavelet paramters
% Assuming wdencmp('gbl',~,~,~,~,~,'h',keepapp);
% with:
%     keepapp = 1

method = 'wavelet';
% https://se.mathworks.com/help/wavelet/ref/wavemngr.html
% wnames = {'sym2', 'bior1.1', 'rbio1.1',...
%     'sym8', 'bior6.8', 'rbio6.8',...
%     'coif1', 'fk4', 'db2', ...
%     'coif5', 'fk22', 'db20' };
% decomp_lvl = [1];
% wtreshold = [0.0005, 0.001, 0.005, 0.01, 0.05];

wnames = {'bior6.8'};
decomp_lvl = 1;
wtreshold = 0.001;


for wn_ii = 1:length(wnames)
    for dl_ii = 1:length(decomp_lvl)
        for cr_ii = 1:length(wtreshold)
            if not(do_wavelet)
                break;
            end
            
            % Set parameter struct
            parameter_struct = struct;
            parameter_struct.method = method;
            parameter_struct.decomp_lvl = decomp_lvl(dl_ii);
            parameter_struct.wname = wnames{wn_ii};
            parameter_struct.wtreshold = wtreshold(cr_ii);
            
            % Store parameter struct
            parameters_cell{end+1} = parameter_struct;
        end
    end
end

%% JPEG2000 Lossy parameters
% JPEG2000 lossy using imwrite

method = 'jpeg2000lossy';
CR_goal = [1.5 2 2.5 3 4 5];

for cr_ii = 1:length(CR_goal)
    if not(do_j2k)
        break;
    end
    
    % Set parameter struct
    parameter_struct = struct;
    parameter_struct.method = method;
    parameter_struct.CR_goal = CR_goal(cr_ii);
    
    % Store parameter struct
    parameters_cell{end+1} = parameter_struct;
end


%% print parameters
if print_check
    for p_set = 1:length(parameters_cell)
        s = parameters_cell(p_set); s = s{1};
        disp(s);
    end
end

%% Save parameters
path = mfilename('fullpath');
rem_str = length(mfilename());
path = path(1:end-rem_str);

save(join([path, 'parameters_cell.mat']), 'parameters_cell')

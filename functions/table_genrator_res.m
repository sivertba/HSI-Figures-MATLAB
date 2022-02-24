clc; close all; clear all;
load('/home/sivertba/OTFP_ENCODING/Data/data_set.mat')
load("results/table.mat");

%% Genearate all table entries

scene_names = ["cuprite", "mofett", "hico"];
bpppbs = [1.0 0.50 0.10 0.05];
methods = ["PCA3D", "PCA2D"];

case_cell_res = {};

%%
counter = 0;
for c_ii = 1:length(case_cell)
    temp_struct = case_cell{c_ii};
    
    if any(contains(methods ,temp_struct.method))
        % do nothing
    else
        continue
    end
    
    data_log = struct;
    
    HSI = data_set.(temp_struct.scene).test(:,:,:);
    scene_name = temp_struct.scene;
    dr_method = temp_struct.method;
    
    data_log.scene = scene_name;
    data_log.dr_method = dr_method;
    
    target_bpppb = temp_struct.bpppb;
    data_log.bpppb = target_bpppb;
    
    target_bpppb = temp_struct.bpppb;
    
    if target_bpppb == 1.0
        target_bpppb_res = 3.0;
    elseif target_bpppb == 0.5
        target_bpppb_res = 3.5;
    elseif target_bpppb == 0.1
        target_bpppb_res = 0.9;
    elseif target_bpppb == 0.05
        target_bpppb_res = 0.95;
    else
        disp("Could not find thing");
        continue        
    end      
    
    best_comps = temp_struct.sam_comps;
    best_j2k_r = temp_struct.sam_j2k_r;
    
    if dr_method == "PCA3D"
        [cube_comp_og,acc_bits_og] = func_j2k_mc_dimred(HSI,best_comps,"PCA",best_j2k_r,scene_name);
    elseif dr_method == "PCA2D"
        [cube_comp_og,acc_bits_og] = func_j2k_dimred(HSI,best_comps,"PCA",best_j2k_r,scene_name);
    elseif dr_method == "OTFP"
        [cube_comp_og,acc_bits_og] = func_j2k_dimred(HSI,best_comps,"OTFP",best_j2k_r,scene_name);
    else
        disp("Could not find thing");
        continue
    end
    
    cube_comp_og = floor(cube_comp_og);
    
    Csize = size(HSI);
    cube_comp_2D = reshape(cube_comp_og,[Csize(1)*Csize(2) Csize(3)]);
    HSI_2D = reshape(HSI,[Csize(1)*Csize(2) Csize(3)]);
    
    cube_comp_og = reshape(cube_comp_2D, Csize);
    
    residuals = HSI_2D - cube_comp_2D;
    
    for res_stdv = 0:1
        
        if res_stdv
            resIndex = WeightedResidualAnalysis(residuals,300,2);
        else
            resIndex = WeightedResidualAnalysis(residuals,1,0);
        end
        sig_res = residuals(resIndex,:);
        
        j2k_r_res = 16/((numel(HSI)/numel(sig_res))*target_bpppb_res);
        
        
        [sig_res_comp, acc_bits_res] = func_j2k_res(sig_res, 1, j2k_r_res, 1);
        
        cube_comp_2D_res = cube_comp_2D;
        cube_comp_2D_res(resIndex,:) = cube_comp_2D_res(resIndex,:) + sig_res_comp;
        
        cube_comp_res = reshape(cube_comp_2D_res, Csize);
        
        % Assuming that the number of elements in HSI image does not exceed 2^32
        bits_res_index_id = size(sig_res_comp,1)+32;
        
        acc_bits_res = acc_bits_res + bits_res_index_id;
        
        bpppb_res = (acc_bits_og+acc_bits_res)/numel(HSI);
        
        data_log.real_bpppb = bpppb_res;
        
        qi_res = QualityIndices(cube_comp_res,HSI);
        
        data_log.res_percent = (numel(HSI)/numel(sig_res));
        data_log.sam = qi_res.sam;
        data_log.snr = qi_res.snr;
        data_log.brisque = qi_res.brisque;
        
        data_log.stdv = res_stdv;
        res_percent = size(sig_res,1)/size(HSI_2D,1);
        
        case_cell_res{end+1} = data_log;
        
    end
    
    save("results/table_res.mat","case_cell_res");
end


%% Plot Case Cell
clc; close all;
load("results/table_res.mat");

data_blocks = struct;

for stdv_ii = [0 1]
    for method_ii = 1:length(methods)
        specific_method = methods(method_ii);
        
        data_block = zeros(length(bpppbs),1 + 3*length(scene_names));
        for case_ii = 1:length(case_cell_res)
            scenario = case_cell_res{case_ii};
            
            ids = find( scene_names == scenario.scene ) -1;
    
            if scenario.dr_method == specific_method && stdv_ii == scenario.stdv
                idx = find(bpppbs == scenario.bpppb);
                if isempty(idx)
                    continue
                end
                
                s_data = data_block(idx,2:end);
                
                try
                    s_data(ids*3+1) = scenario.snr;
                    s_data(ids*3+2) = scenario.sam;
                    s_data(ids*3+3) = scenario.brisque;
                catch
                    % do nothing
                end
                
                data_block(idx, :) = [scenario.bpppb, s_data];
            end
        end

        data_blocks.(specific_method + num2str(stdv_ii)) = data_block;
    end
end

save("results/data_blocks_res.mat","data_blocks");

%% Plot per bpppb
clc; close all;
load("results/data_blocks_res.mat");

fileID = fopen("results/bppps_res.csv",'w');

fprintf(fileID, strjoin(["bpppb" "MET" repmat(scene_names(1),1,3) repmat(scene_names(2),1,3) repmat(scene_names(3),1,3)],',')+"\n");
fprintf(fileID, strjoin(["     " "     " repmat(["SNR" "SAM" "BQ"],1,3)], ','));
fprintf(fileID, '\n');

db_fieldnames = fieldnames(data_blocks);
for rate_ii = 1:length(bpppbs)    
    for met_ii = 1:length(db_fieldnames)
        if db_fieldnames{met_ii} == "NoSpat"
            continue;
        end
        block = data_blocks.(db_fieldnames{met_ii});
        idb = find( block(:,1) == bpppbs(rate_ii));
        
        fprintf(fileID, "%.2f,%s", bpppbs(rate_ii), db_fieldnames{met_ii});
        
        for z_ii = 2:size(block,2)
            fprintf(fileID, ',%.2f', block(idb, z_ii));
        end
        fprintf(fileID, "\n");
    end
end

fclose(fileID);

clc; close all; clear all;
load('/home/sivertba/OTFP_ENCODING/Data/data_set.mat')

%% Genearate all table entries
case_cell = {};

scene_names = ["cuprite", "mofett", "hico"];
bpppbs = [4.0 1.0 0.50 0.10 0.05];
methods = ["None3D", "PCA3D", "None2D", "PCA2D", "OTFP","NoSpat"];

for s_ii = 1:length(scene_names)
    for b_ii = 1:length(bpppbs)
        for m_ii = 1:length(methods)
            case_struct = struct;
            case_struct.scene = scene_names(s_ii);
            case_struct.bpppb = bpppbs(b_ii);
            case_struct.method = methods(m_ii);
            
            case_cell{end+1} = case_struct;
        end
    end
end

%%
for c_ii = 1:length(case_cell)
    temp_struct = case_cell{c_ii};
    disp(join(["Case" c_ii " of " length(case_cell)]));
    disp(temp_struct);
    
    HSI = data_set.(temp_struct.scene).test(:,:,:);
    scene_name = temp_struct.scene;
    dr_method = temp_struct.method;
    
    target_bpppb = temp_struct.bpppb;
    
    if contains(dr_method ,"None",'IgnoreCase',true)
        comps = size(HSI,3);
        j2k_r = 16/target_bpppb;
        
    else
        comps = 1:size(HSI,3);
        j2k_r = (16/target_bpppb)*comps/size(HSI,3);
        
        comps = comps(j2k_r > 3);
        j2k_r = j2k_r(j2k_r > 3);
    end
    
    
    best_sam = 200;
    best_snr = 0;
    for it = 1:length(j2k_r)
        if dr_method == "NoSpat"
            [cube_comp,acc_bits] = func_j2k_dimred_no_spat(HSI,comps(it),"PCA",j2k_r(it),scene_name);
        elseif dr_method == "OTFP3D"
            [cube_comp,acc_bits] = func_j2k_mc_dimred(HSI,comps(it),"OTFP",j2k_r(it), scene_name);
        elseif dr_method == "PCA3D"
            [cube_comp,acc_bits] = func_j2k_mc_dimred(HSI,comps(it),"PCA",j2k_r(it), scene_name);
        elseif dr_method == "None3D"
            [cube_comp,acc_bits] = func_j2k_mc_dimred(HSI,comps(it),"None",j2k_r(it), scene_name);
        elseif dr_method == "None2D"
            [cube_comp,acc_bits] = func_j2k_dimred(HSI,comps(it),"None",j2k_r(it),scene_name);
        elseif dr_method == "PCA2D"
            [cube_comp,acc_bits] = func_j2k_dimred(HSI,comps(it),"PCA",j2k_r(it),scene_name);
        elseif dr_method == "OTFP"
            [cube_comp,acc_bits] = func_j2k_dimred(HSI,comps(it),"OTFP",j2k_r(it),scene_name);
        end
        
        qi = QualityIndices(cube_comp,HSI);
        
        if qi.sam < best_sam
            best_sam = qi.sam;
            case_cell{c_ii}.sam_sam = qi.sam;
            case_cell{c_ii}.sam_bq = qi.brisque;
            case_cell{c_ii}.sam_snr = qi.snr;
            case_cell{c_ii}.sam_comps = comps(it);
            case_cell{c_ii}.sam_j2k_r = j2k_r(it);
            case_cell{c_ii}.sam_bits = acc_bits;
            case_cell{c_ii}.sam_real_bpppb= acc_bits/numel(HSI);
        end
        
        if qi.snr > best_snr
            best_snr = qi.snr;
            case_cell{c_ii}.snr_sam = qi.sam;
            case_cell{c_ii}.snr_bq = qi.brisque;
            case_cell{c_ii}.snr_snr = qi.snr;
            case_cell{c_ii}.snr_comps = comps(it);
            case_cell{c_ii}.snr_j2k_r = j2k_r(it);
            case_cell{c_ii}.snr_bits = acc_bits;
            case_cell{c_ii}.snr_real_bpppb= acc_bits/numel(HSI);
        end
    end
    save("results/table.mat","case_cell","c_ii");
end


%% Plot Case Cell
clc; close all;
load("results/table.mat");

data_blocks = struct;

min_j2kr = [100 100 100];
min_comps = [200 200 200];
bands = [188 188 87];

for method_ii = 1:length(methods)
    specific_method = methods(method_ii);
    
    data_block = zeros(length(bpppbs),1 + 3*length(scene_names));
    for case_ii = 1:length(case_cell)
        scenario = case_cell{case_ii};
        
        ids = find( scene_names == scenario.scene ) -1;
        
        if min_j2kr(ids+1) > scenario.sam_j2k_r
            min_j2kr(ids+1) = scenario.sam_j2k_r;
        end
        
        re_comps = scenario.sam_j2k_r*scenario.bpppb*bands(ids+1)/16;
        if min_comps(ids+1) > re_comps
            min_comps(ids+1) = re_comps;
        end
        
        if scenario.method == specific_method
            idx = find(bpppbs == scenario.bpppb);
            if isempty(idx)
                continue
            end
            
            s_data = data_block(idx,2:end);
            
            try
                s_data(ids*3+1) = scenario.sam_snr;
                s_data(ids*3+2) = scenario.sam_sam;
                s_data(ids*3+3) = scenario.sam_bq;
            catch
                % do nothing
            end
            
            data_block(idx, :) = [scenario.bpppb, s_data];
        end
    end
    data_blocks.(strrep(specific_method,'-','')) = data_block;
end

save("results/data_blocks.mat","data_blocks");

%% Plot per bpppb
clc; close all;
load("results/table.mat");

fileID = fopen("results/bppps_plain.csv",'w');

fprintf(fileID, strjoin(["bpppb" "MET" repmat(scene_names(1),1,3) repmat(scene_names(2),1,3) repmat(scene_names(3),1,3)],',')+"\n");
fprintf(fileID, strjoin(["     " "     " repmat(["SNR" "SAM" "BQ"],1,3)], ','));

db_fieldnames = fieldnames(data_blocks);
for rate_ii = 1:length(bpppbs)
    fprintf(fileID, '%.2f\n', bpppbs(rate_ii));
    disp("");
    
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

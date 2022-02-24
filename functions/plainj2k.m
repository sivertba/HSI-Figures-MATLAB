clc; close all; clear all;
load('/home/sivertba/OTFP_ENCODING/Data/data_set.mat')
%%
with_residuals = 1;

HSI = data_set.mofett.test(:,:,:);
scene_name = "mofett";

comps = 176;
dr_method = "PCA";

target_bpppb = 1.0;
res_part = 0.1;

j2k_r = (16/target_bpppb)*comps/size(HSI,3);

[cube_comp,acc_bits] = func_j2k_dimred(HSI,comps,dr_method,j2k_r,scene_name);
% [cube_comp,acc_bits] = func_j2k_dimred_no_spat(HSI,comps,dr_method,j2k_r,scene_name);
% [cube_comp,acc_bits] = func_j2k_mc_dimred(HSI,comps,dr_method,j2k_r, scene_name);

if with_residuals
    target_bpppb_og = (1-res_part)*target_bpppb;
    target_bpppb_res = res_part*target_bpppb;

    j2k_r_og = (16/target_bpppb_og)*comps/size(HSI,3);
    
    [cube_comp_og,acc_bits_og] = func_j2k_dimred(HSI,comps,dr_method,j2k_r_og,scene_name);
    
    Csize = size(HSI);
    cube_comp_2D = reshape(cube_comp_og,[Csize(1)*Csize(2) Csize(3)]);
    HSI_2D = reshape(HSI,[Csize(1)*Csize(2) Csize(3)]);
        
    cube_comp_og = reshape(cube_comp_2D, Csize);

    residuals = HSI_2D - cube_comp_2D;
    
    resIndex = WeightedResidualAnalysis(residuals,200,3);
    sig_res = residuals(resIndex,:);
    
    j2k_r_res = 16/((numel(HSI)/numel(sig_res))*target_bpppb_res);
    [sig_res_comp, acc_bits_res] = func_j2k_res(sig_res, 1, j2k_r_res);
%     [sig_res_comp, acc_bits_res] = func_j2k_res(sig_res, 1, 1);
    
    cube_comp_2D_res = cube_comp_2D;
    cube_comp_2D_res(resIndex,:) = cube_comp_2D_res(resIndex,:) + sig_res_comp;
    
    cube_comp_res = reshape(cube_comp_2D_res, Csize);
    
    % Assuming that the number of elements in HSI image does not exceed 2^32
    bits_res_index_id = size(sig_res_comp,1)+32; 
    
    acc_bits_res = acc_bits_res + bits_res_index_id;
    
end

bpppb = acc_bits/numel(HSI);
disp(join(["bpppb:" bpppb]));

qi = QualityIndices(cube_comp,HSI);
disp(qi);

if with_residuals
    bpppb_res = (acc_bits_og+acc_bits_res)/numel(HSI);
    disp(join(["bpppb with residuals:" bpppb_res]));
    
    qi_res = QualityIndices(cube_comp_res,HSI);
    disp(qi_res);
end

hsiplay = @(M,Msad) implay([M/max(M(:)) Msad/max(M(:)) 100*(M-Msad)/max(M(:))]);
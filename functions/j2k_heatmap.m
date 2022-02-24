clc; close all; clear all;
load('/home/sivertba/OTFP_ENCODING/Data/data_set.mat')
%%
HSI = data_set.cuprite.test(:,:,:);

dr_method = "PCA";

CR_map = zeros(30,30);
snr_map = zeros(30,30);


for comps = 1:30
    for j2k_rs = 1:30
        
        [cube_comp,acc_bits] = func_j2k_dimred(HSI,comps,dr_method,j2k_rs);
        
        
        CR_map(comps,j2k_rs) = numel(HSI)*16 / acc_bits;
        
        [snr_val, ~] = plainSNR(cube_comp,HSI);
        snr_map(comps,j2k_rs) = snr_val.snr_2;
        
    end
end




clc; close all; clear all;
load('/home/sivertba/OTFP_ENCODING/Data/data_set.mat')


%% converting

methods = ["PCA2D"];
T_s1 = table();
T_s2 = table();

target_bpppbs = 4;

d_set = ["cuprite" "mofett"];

for d_ii = 1:length(d_set)
    
    target_bpppb = target_bpppbs;
    
    HSI = data_set.(d_set(d_ii)).test;
    Osize = size(HSI);
    X_test= reshape(HSI,[Osize(1)*Osize(2) Osize(3)]);
    
    dr_method = methods;
    
    comp = 80;
    
    HSI_train = data_set.(d_set(d_ii)).tl;
    Tsize = size(HSI_train);
    X_train= reshape(HSI_train,[Tsize(1)*Tsize(2) Tsize(3)]);
    
    loadings_pca = pca(X_train);
    loadings = loadings_pca(:,:)';
    
    scores_res = loadings(1:comp,:)*X_test';
    
    tic;
    scores = loadings(1:comp,:)*X_test';
    spatialized_cube = spatialize_scores(scores', size(HSI));
    
    
    bits = 2^16 - 1;
    t_min=min(spatialized_cube(:));
    t_max=max(spatialized_cube(:));
    t_new = round(rescale(spatialized_cube,0,bits));
    
    fname = 'll_1.j2k';
    imwrite(uint16(t_new), fname, 'mode','lossless');
    time_s1 = toc;
    
    acc_bits = 0;
    s=dir(fname);
    
    acc_bits = acc_bits + s.bytes*8;
    
    spatialized_cube_comp = rescale(imread(fname),t_min,t_max);
    scores_comp = transpose(despatialize_scores(spatialized_cube_comp, size(HSI), size(scores(1:comp,:)')));
    X_comp = scores_comp(1:comp,:)' * loadings(1:comp,:);
    
    cube_comp = reshape(X_comp,[Osize(1) Osize(2) Osize(3)]);
    
    qi_reg = QualityIndices(cube_comp,HSI);
    
    
    tic;
    
    residuals = HSI-round(cube_comp);
    
    floor_val = min(residuals(:));
    res_rise = residuals - floor_val;
    
    fname_res = 'resll.j2k';
    
    imwrite(res_rise, fname_res, 'mode','lossless');
    time_s2 = toc;
    
    acc_bits_res = 0;
    s=dir(fname_res);
    acc_bits_res = acc_bits_res + s.bytes*8;
    
    read_residuals = double(imread(fname));
    res_out = read_residuals+floor_val;
    
    cube_comp_res = read_residuals+cube_comp;
    
    bpppb_fin = (acc_bits+acc_bits_res)/numel(HSI);
    
    qi_res = QualityIndices(cube_comp_res,HSI);
    
    C = {dr_method,...
        scene,...
        bpppb_fin,...
        comp,...
        qi_reg.snr,...
        qi_reg.sam,...
        qi_reg.brisque,...
        time_s1};
    T_new = cell2table(C,...
        'VariableNames',...
        {'method' 'scene' 'bpppb_fin' 'comps' 'snr' 'sam' 'bq' 'time'});
    T_s2 = [T_s2; T_new];
end

%% Save results
writetable(T_s2,"figs/lossless.csv");

%%

h_data = data_set.hico.test;

h_data_ish = round(h_data*50);

imwrite(uint16(h_data_ish), "test.j2k", 'mode','lossless');
fout = imread("test.j2k");
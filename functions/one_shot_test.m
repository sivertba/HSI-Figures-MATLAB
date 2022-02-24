clc; close all; clear all;
load('/home/sivertba/OTFP_ENCODING/Data/data_set.mat')

scene = "Maine";
scene_name = scene;


%% converting

methods = ["PCA2D"];
T_s1 = table();
T_s2 = table();

% comps = 18:18:180;
comps = 87;
% comps = 45:45:180; % test

target_bpppbs = 2;
target_bpppb_res_li = 2;

for t_bppps = 1:length(target_bpppbs)
    for comp_ii = 1:length(comps)
        
        target_bpppb = target_bpppbs(t_bppps);
        
%         HSI = data_set.hico.test;
        HSI = out;
        Osize = size(HSI);
        X_test= reshape(HSI,[Osize(1)*Osize(2) Osize(3)]);
        
        dr_method = methods;
        
        comp = comps(comp_ii);
        best_j2k_r = (16/target_bpppb)*comp/Osize(3);
        
        
        HSI_train = data_set.hico.tl;
        HSI_train = out;
        Tsize = size(HSI_train);
        X_train= reshape(HSI_train,[Tsize(1)*Tsize(2) Tsize(3)]);
        
        mu = mean(X_train,1);
        
        X_train = X_train - mu;
        X_test = X_test - mu;
        
        loadings_pca = pca(X_train, 'Algorithm', 'eig');
        loadings = loadings_pca(:,:)';
        
        scores_res = loadings(1:comp,:)*X_test';
        residuals = X_test - scores_res(1:comp,:)' * loadings(1:comp,:);
        
        tic;
        scores = loadings(1:comp,:)*X_test';
        spatialized_cube = spatialize_scores(scores', size(HSI));
        
        system('rm j2k_test/*');
        
        bits = 2^16 - 1;
        t_min=min(spatialized_cube(:));
        t_max=max(spatialized_cube(:));
        t_new = round(rescale(spatialized_cube,0,bits));
        
        imwrite(uint16(t_new), 'j2k_test/png_1.png');
        system(join(['opj_compress -ImgDir j2k_test/ -OutFor J2K -r ' num2str(best_j2k_r)]));
        time_s1 = toc;
        
        system('for i in j2k_test/*.J2K; do opj_decompress -i $i -o j2k_test/$(basename "$i" .J2K)_rev.png ; done;');
        
        acc_bits = 0;
        s=dir('j2k_test/png_1.J2K');
        
        acc_bits = acc_bits + s.bytes*8;
        
        fname = 'j2k_test/png_1_rev.png';
        
        spatialized_cube_comp = rescale(imread(fname),t_min,t_max);
        scores_comp = transpose(despatialize_scores(spatialized_cube_comp, size(HSI), size(scores(1:comp,:)')));
        X_comp = scores_comp(1:comp,:)' * loadings(1:comp,:) + mu;
        
        cube_comp = reshape(X_comp,[Osize(1) Osize(2) Osize(3)]);
        
        qi_reg = QualityIndices(cube_comp,HSI);
        
        C = {dr_method,...
            scene,...
            target_bpppb,...
            comp,...
            qi_reg.snr,...
            qi_reg.sam,...
            qi_reg.brisque,...
            time_s1};
        T_new = cell2table(C,...
            'VariableNames',...
            {'method' 'scene' 'bpppb' 'comps' 'snr' 'sam' 'bq' 'time'});
        
        disp(T_new);
        T_s1 = [T_s1; T_new];
        
        target_bpppb_res = target_bpppb_res_li(t_bppps);
        
        if target_bpppb_res
            tic;
            
%             resIndex = ones(size(residuals, 1), 1);
                                    
            mu = mean(residuals);
            sigma = std(residuals);
            nStdv = 1;
            
            resIndex = zeros(size(residuals, 1), 1);
            for bandInd = 1:size(residuals, 2)
                resLB = mu(bandInd) - sigma(bandInd) * nStdv;
                resUB = mu(bandInd) + sigma(bandInd) * nStdv;
                resIndex = resIndex | ...
                    residuals(:, bandInd) < resLB | ...
                    residuals(:, bandInd) > resUB;
            end
            
            sig_res = residuals(resIndex,:);
            non_sig_res = resIndex - ones(size(residuals, 1), 1);

        end
    end
end
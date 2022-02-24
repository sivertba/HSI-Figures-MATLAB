%% fakk
clc; close all; clear all;
load('/home/sivertba/OTFP_ENCODING/Data/data_set.mat')
hsiplay = @(M,Msad) implay([M/max(M(:)); Msad/max(M(:)); 10*(M-Msad)/max(M(:))]);

%%
HSI = data_set.mofett.test;
Osize = size(HSI);
X_test= reshape(HSI,[Osize(1)*Osize(2) Osize(3)]);

comp = 185;        

% loadings_pca = pca(X_test);
% loadings = loadings_pca(:,:)';
%     

loadings = pca(X_test);
loadings = loadings';
scores = loadings(:,1:comp)'*X_test';

residuals = X_test - scores(1:comp,:)' * loadings(1:comp,:);

times = zeros(100,2);
for ii = 1:100
    tic;
    resIndex1 = WeightedResidualAnalysis(residuals,100,2);
    idle_time = toc;
    
    tic;

    mu = mean(residuals);
    sigma = std(residuals);
    nStdv = 2;
    
    identifiedRows = zeros(size(residuals, 1), 1);
    for bandInd = 1:size(residuals, 2)
        resLB = mu(bandInd) - sigma(bandInd) * nStdv;
        resUB = mu(bandInd) + sigma(bandInd) * nStdv;
        identifiedRows = identifiedRows | ...
            residuals(:, bandInd) < resLB | ...
            residuals(:, bandInd) > resUB;
    end
    alt_time = toc;

    times(ii,:) = [idle_time alt_time];
end


X_comp = scores'* loadings(:,1:comp)';
HSI_new= reshape(X_comp,[Osize(1) Osize(2) Osize(3)]);



hsiplay(HSI,HSI_new)



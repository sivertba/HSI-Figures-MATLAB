clc; close all; clear;
%% Test for score spatialization

% Load data
load('hico_data.mat')

% choose which scans to use from datacube (nPixels x nScans x nChannels)
firstScan = 1;
nScans = 320;

rad = test(:,firstScan:firstScan+nScans - 1,1:87);

cube_size = size(rad);
X = reshape(rad,[cube_size(1)*cube_size(2) cube_size(3)]);
Nindex=1:(cube_size(1)*cube_size(2));

disp('Data and parameters loaded ... ')

% Data parameters
% scoreThresh = zeros(3,20);
% totalWeight = ones(nScans,87);
% wpi = FindRefPixel(X);
% I0 = mean(X(wpi,:));
% totalWeight = totalWeight./I0;
% 
% center=mean(X./I0);
% X = PreprocessRawData(X, center, totalWeight, cube_size(2));

% Make scores
[~,loadings] = dimred_func(X,'PCA');
[scores, ~] = ProjectData(X, loadings(1:3,:));

% Spatialize scores and show results
spatialized_scores = spatialize_scores(scores, cube_size);

fig1 = figure;

imshow(create_RGB(test,38,23,7,0.08,0.5));
exportgraphics(fig1,'figs/hico_rgb.png')

fig2 = figure;
imshow(spatialized_scores/max(spatialized_scores(:)));
exportgraphics(fig2,'figs/hico_3scores.png')
%%
fig3 = figure;

[coeff,~,latent] = pca(X);
s_coeff =[coeff(:,1)*latent(1) coeff(:,2)*latent(2) coeff(:,3)*latent(3)];
wl_hico = [ 404.080 409.808 415.536 421.264 426.992 432.720 438.448	 444.176 449.904 455.632 461.360 467.088 472.816 478.544 484.272 490.000 495.728 501.456 507.184 512.912 518.640 524.368 530.096	 535.824 541.552 547.280 553.008 558.736 564.464 570.192 575.920 581.648 587.376 593.104 598.832 604.560 610.288 616.016 621.744 627.472 633.200 638.928 644.656 650.384 656.112 661.840 667.568 673.296 679.024 684.752 690.480 696.208 701.936 707.664 713.392 719.120 724.848 730.576 736.304 742.032 747.760 753.488 759.216 764.944 770.672 776.400 782.128 787.856 793.584 799.312 805.040 810.768 816.496 822.224 827.952 833.680 839.408 845.136 850.864 856.592 862.320 868.048 873.776 879.504 885.232 890.960 896.688];
plot(wl_hico, s_coeff);
grid on;

legend("Component 1","Component 2","Component 3",'Location',"best");
ylabel("Magnitude");
xlabel("Wavelength (nanometer)");
title("First Principal Components of HICO");
set(fig3,'position',[104 600 833 310])


exportgraphics(fig3,'figs/hico_loadings.png')


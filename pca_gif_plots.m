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
scoreThresh = zeros(3,20);
totalWeight = ones(nScans,87);
wpi = FindRefPixel(X);
I0 = mean(X(wpi,:));
totalWeight = totalWeight./I0;

center=mean(X./I0);
X = PreprocessRawData(X, center, totalWeight, cube_size(2));

% Make scores
[~,loadings] = dimred_func(X,'PCA');
[scores, ~] = ProjectData(X, loadings(:,:));

% Spatialize scores and show results
spatialized_scores = spatialize_scores(scores, cube_size);


%% Get Spatial PNGs
max_val = max(spatialized_scores(:));
img_li =1:size(spatialized_scores,1):size(spatialized_scores,2);
img_li_end = img_li +size(spatialized_scores,1)-1;

for ii = 1:size(test,3)
    fig = figure;
    
    frame = spatialized_scores(:,img_li(ii):img_li_end(ii));
    
    subplot(121);
    imshow(frame/max_val);
    title("Max Normalized")
    
    subplot(122);
    frame = spatialized_scores(:,img_li(ii):img_li_end(ii));
    imshow(frame/max(frame(:)));
    title("Comp Normalized")
    
    box_string = sprintf("c:    %s/%s",num2str(ii,'%02.f'),num2str(size(test,3)));
    a = annotation('textbox',[.45 .33 .4 .5],'String',box_string,'FitBoxToText','on');
    a.Color = 'w';
    a.BackgroundColor = 'k';
    
%     fname = sprintf('figs/pca_gif/pca_gif%s.png',num2str(ii,'%02.f'));
%     exportgraphics(fig,fname)
    plot_darkmode;
    break;
    
end

%%
close all;


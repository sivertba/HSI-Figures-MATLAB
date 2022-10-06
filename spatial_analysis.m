%% spatial analysis
clc; clear; close all;
mkdir("output/");
mkdir("output/spatial/")

%% Load data
load("h1_data.mat")
scene_name = "kuwait";

%% define spatial ROI
xstart  = 140;
xend    = 240;
ystart  = 450;
yend    = 650;
spatial_roi_img = cube(xstart:xend,ystart:yend,:);

%% Plot and store images
for bb = 1:size(cube,3)
    imagesc(spatial_roi_img(:,:,bb));

    box_string = sprintf("%s nm",num2str(h1_wl(bb),'%02.f'));
    a = annotation('textbox',[.8 .1 .15 .1],'String',box_string,'FitBoxToText','on');
    a.Color = 'w';
    a.BackgroundColor = 'k';
    
    fname = sprintf('output/spatial/%s_%snm.png',scene_name, num2str(h1_wl(bb),'%02.f'));
    exportgraphics(gcf,fname,"Resolution",400);
end

%% Plot Line 
I = spatial_roi_img(:,:,54);
imagesc(I);
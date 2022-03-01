clc; close all; clear;

% Load data
% data = "hico"
data = "h1";

if data == "hico"
    load('hico_data.mat')

    wl_hico = [ 404.080 409.808 415.536 421.264 426.992 432.720 438.448	 444.176 449.904 455.632 461.360 467.088 472.816 478.544 484.272 490.000 495.728 501.456 507.184 512.912 518.640 524.368 530.096	 535.824 541.552 547.280 553.008 558.736 564.464 570.192 575.920 581.648 587.376 593.104 598.832 604.560 610.288 616.016 621.744 627.472 633.200 638.928 644.656 650.384 656.112 661.840 667.568 673.296 679.024 684.752 690.480 696.208 701.936 707.664 713.392 719.120 724.848 730.576 736.304 742.032 747.760 753.488 759.216 764.944 770.672 776.400 782.128 787.856 793.584 799.312 805.040 810.768 816.496 822.224 827.952 833.680 839.408 845.136 850.864 856.592 862.320 868.048 873.776 879.504 885.232 890.960 896.688];
    wl = round(wl_hico);

else
    load("h1_data.mat");
    test = h1_data;
    wl = round(h1_wl);
end

%% Get Spatial PNGs
for ii = 1:size(test,3)
    fig = figure;
    imshow(test(:,:,ii)/max(test(:)));
    box_string = sprintf("nm:  %s\nf:      %s/%s",...
        num2str(wl(ii)),num2str(ii,'%02.f'),num2str(size(test,3)));
    a = annotation('textbox',[.2 .425 .4 .5],'String',box_string,'FitBoxToText','on');
    a.Color = 'w';
    a.BackgroundColor = 'k';
    
    
    fname = sprintf('figs/spatial_gif/spatial%s_gif%s.png',data,num2str(ii,'%02.f'));
%     fig.Position = [1636 674 700 600];
    exportgraphics(fig,fname)
%     pause(.5)
%     break;
end

%%
close all;
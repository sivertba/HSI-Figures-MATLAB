clc; close all; clear;
%% Test for score spatialization

% Load data
load('hico_data.mat')

wl_hico = [ 404.080 409.808 415.536 421.264 426.992 432.720 438.448	 444.176 449.904 455.632 461.360 467.088 472.816 478.544 484.272 490.000 495.728 501.456 507.184 512.912 518.640 524.368 530.096	 535.824 541.552 547.280 553.008 558.736 564.464 570.192 575.920 581.648 587.376 593.104 598.832 604.560 610.288 616.016 621.744 627.472 633.200 638.928 644.656 650.384 656.112 661.840 667.568 673.296 679.024 684.752 690.480 696.208 701.936 707.664 713.392 719.120 724.848 730.576 736.304 742.032 747.760 753.488 759.216 764.944 770.672 776.400 782.128 787.856 793.584 799.312 805.040 810.768 816.496 822.224 827.952 833.680 839.408 845.136 850.864 856.592 862.320 868.048 873.776 879.504 885.232 890.960 896.688];
wl_hico = round(wl_hico);


cube_size = size(test);
X = reshape(test,[cube_size(1)*cube_size(2) cube_size(3)]);

% Select orginal spectra
random_id = 1000;
orginal_spectra = X(random_id,:);


% Make scores
[~,loadings] = dimred_func(X,'PCA');
%%
% orginal_spectra = X(random_id,:);
% new_spectra = orginal_spectra*pinv(loadings(1:ii,:))*loadings(1:ii,:);
% plot(wl_hico, orginal_spectra, wl_hico,new_spectra);

for ii = 1:size(test,3)
    fig = figure;

    orginal_spectra = X(random_id,:);
    new_spectra = orginal_spectra*pinv(loadings(1:ii,:))*loadings(1:ii,:);
    plot(wl_hico, orginal_spectra,'r', wl_hico,new_spectra,'m');
    
    box_string = sprintf("Components:    %s/%s",num2str(ii,'%02.f'),num2str(size(test,3)));
    a = annotation('textbox',[.15 .8 .1 .1],'String',box_string,'FitBoxToText','on');
    a.Color = '#edb99c';
    a.BackgroundColor = 'k';

    xlabel("Wavelength (nm)")
    ylabel("Radiance (W/m2/um/sr)")

    set(gca,'Color','k')
    set(gca,'XColor','#edb99c')
    set(gca,'YColor','#edb99c')

    legend(["\color{red}Original Spectra", "\color{magenta}Resotred Spectra"]);


    axis([ 400 900 -5 45])
    
    fname = sprintf('figs/spectra_gif/spectra_gif%s.png',num2str(ii,'%02.f'));
    exportgraphics(fig,fname)
%     break;
    
end

%%
close all;


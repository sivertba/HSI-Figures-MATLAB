function [QI_COMP, QI_LOAD] = plot_results_func(rad_ORIG,rad_LOAD,...
    rad_COMP, analysis_results,compression_parameters,...
    param_set, scene, dr_method)
%% Plot results

% unpack analysis_results
compressedSize = analysis_results.compressedSize;
compressedSizewoRes = analysis_results.compressedSizewoRes;
EV = analysis_results.EV;
EVwoRes = analysis_results.EVwoRes;

compressedSize_comp = analysis_results.compressedSize_comp;
compressedSizewoRes_comp = analysis_results.compressedSizewoRes_comp;
EV_comp = analysis_results.EV_comp;
EVwoRes_comp = analysis_results.EVwoRes_comp;

% QI Metrics
QI_ORIG = QualityIndices(rad_ORIG,rad_ORIG);
QI_LOAD = QualityIndices(rad_LOAD,rad_ORIG);
QI_COMP = QualityIndices(rad_COMP,rad_ORIG);

fig_handle = figure;

cp_method = compression_parameters.method;

if cp_method == "fft"
    fft_tresh = compression_parameters.tresh;
    sgtitle(scene+' - '+dr_method+' - FFT, w/ tresh of '+ string(fft_tresh));

elseif cp_method == "wavelet"
    decomp_lvl = compression_parameters.decomp_lvl;
    wname = compression_parameters.wname;
    wtreshold = compression_parameters.wtreshold;
    
    sgtitle(join([scene " - " dr_method " - Wavelet, with:" "decomp lvl: "...
        decomp_lvl ", wname: " wname ", wtreshold: " wtreshold]));

elseif cp_method == "jpeg2000lossy"
    cr_goal = compression_parameters.CR_goal;
    sgtitle(scene+' - '+dr_method+' - JPEG2000, w/ CR goal of '+ string(cr_goal));

else
    error([cp_method " is not a known method for plotting"]);
end

hold on;

% original
subplot(3,3,1);
imshow(create_RGB(rad_ORIG,38,23,7,0.08,0.5))
title("Full Image");

subplot(3,3,3+1);
imshow(create_RGB(rad_ORIG(250:1:300,250:1:300,:),38,23,7,0.08,0.5))
title("Zoomed Image");

ax = subplot(3,3,6+1);
str_ORIG = sprintf("Original Image \nNo Filter\n\n"+...
    "Compression ratio w/o residuals = 1 \n"+...
    "Mean bandwidth PSNR: = %.4f\n"+...
    "Cross Correlation: = %.4f\n"+...
    "Spectral Angle Map: = %.4f\n"+...
    "Normalized RMSE: = %.4f\n",...
    QI_ORIG.psnr,QI_ORIG.cc,QI_ORIG.sam, QI_ORIG.nrmse);
text(0.5,0.5,str_ORIG,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle');
title("Comparison Summary");
set ( ax, 'visible', 'off')

% LOAD
subplot(3,3,2);
imshow(create_RGB(rad_LOAD,38,23,7,0.08,0.5))
title("Full Image");

subplot(3,3,3+2);
imshow(create_RGB(rad_LOAD(250:1:300,250:1:300,:),38,23,7,0.08,0.5))
title("Zoomed Image");

ax = subplot(3,3,6+2);
str_LOAD = sprintf("Restored Image \nUsing 20 Loadings\n\n"+...
    "Compression ratio w/o residuals = %.4f \n"+...
    "Mean bandwidth PSNR: = %.4f\n"+...
    "Cross Correlation: = %.4f\n"+...
    "Spectral Angle Map: = %.4f\n"+...
    "Normalized RMSE: = %.4f\n",...
    1/compressedSizewoRes, ...
    QI_LOAD.psnr,QI_LOAD.cc,QI_LOAD.sam, QI_LOAD.nrmse);
text(0.5,0.5,str_LOAD,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle');
title("Comparison Summary");
set ( ax, 'visible', 'off')

% COMP
subplot(3,3,3);
imshow(create_RGB(rad_COMP,38,23,7,0.08,0.5))
title("Full Image");

subplot(3,3,3+3);
imshow(create_RGB(rad_COMP(250:1:300,250:1:300,:),38,23,7,0.08,0.5))
title("Zoomed Image");

ax = subplot(3,3,6+3);
str_COMP  = sprintf("Restored Image \nAfter compression\n\n"+...
    "Compression ratio w/o residuals = %.4f \n"+...
    "Mean bandwidth PSNR: = %.4f\n"+...
    "Cross Correlation: = %.4f\n"+...
    "Spectral Angle Map: = %.4f\n"+...
    "Normalized RMSE: = %.4f\n",...
    1/compressedSizewoRes_comp,...
    QI_COMP.psnr,QI_COMP.cc,QI_COMP.sam, QI_COMP.nrmse);
text(0.5,0.5,str_COMP,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle');
title("Comparison Summary");
set ( ax, 'visible', 'off')

set(fig_handle,'position',[200, 300, 1300, 1000])

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir results;

saveas(fig_handle,...
    ['results/vd_' cp_method '_id' num2str(param_set) '_' date '.png'],'png');

close all;
end


%%  Wavelet compressed scores and percentage of thresholded coefficients 
%   i.e. coefficients that are not equal to 0
function [compressed_scores, thresholded_coefficients, transform] = wavelet_compression(original_scores, parameters)
% % Input:  original_scores, A numeric matrix representing the scores.
% % Input:  parameters, A struct structure with parameters for wavelets
% %             - method, string saying wavelet.
% %             - decomp_lvl, Decomposition level specified as a 
% %               positive integer. wavedec2 does not enforce a maximum 
% %               level restriction. Use wmaxlev to determine the maximum 
% %               decomposition level possible of the matrix X using the 
% %               wavelet wname. The maximum level is the last level for 
% %               which at least one coefficient is correct.
% %             - wname, Analyzing wavelet, specified as a character 
% %               vector or string scalar.
% %             - wtreshold, Threshold to apply to the wavelet coefficients, 
% %               specified as a scalar, real-valued vector, or 
% %               real-valued matrix. For the case 'gbl', 
% %               wtreshold is a scalar.
% % 
% % Output: compressed_scores, A numeric matrix representing the scores of 
% %         after compression with wavelets.
% % Output: thresholded_coefficients, A scalar representing the Compression 
% %         score, the percentage of thresholded coefficients that are not 0.
% % Output: transform, A struct containing the variables needed for the 
% %         transformation between wavelet space and orgiinal


assert(parameters.method == "wavelet", "Wrong set of parameters");
% assert(wmaxlev(original_scores,parameters.wname) > parameters.decomp_lvl...
%     ,"Decomposition level is too high for the number of scores");

% decomp_lvl = parameters.decomp_lvl;

wname = parameters.wname;
decomp_lvl = wmaxlev(size(original_scores),wname);
wtreshold = parameters.wtreshold;

transform = struct;

[c,l] = wavedec2(original_scores,decomp_lvl,wname);



keepapp = 0;
[compressed_scores,c_new,~,perf0,~] = ...
    wdencmp('gbl',c,l,wname,decomp_lvl,wtreshold,'h',keepapp);

transform.wl_decomp = c_new;
transform.wl_bookkeeping = l;
transform.wl_lvl = decomp_lvl;
transform.wl_wname= wname;

thresholded_coefficients = (100-perf0)/100;

end
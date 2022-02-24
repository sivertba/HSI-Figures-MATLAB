function [post_residuals_2d,acc_bits] = func_j2k_res(residuals_2d, decorrelation, j2k_r, do_rescale)
% %
% % Input residuals_2d: image cube to be compressed (spatial x spatial x spectral)
% %
% % Input decorrelation: wheter or not to use spectral decorrelation
% %
% % Input j2k_r: compression ratio set by JPEG2000 stage
% %
% % Output res_comp: Restored residuals after decompression.
% %
% % Output acc_bits: Number of bits after compression.

%%%%

if decorrelation
    [scores, loadings] = dimred_func(residuals_2d,"PCA");
    res_img = scores;
else
    res_img = residuals_2d;
end

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir j2k_res;

if do_rescale
    t_min=min(res_img(:));
    t_max=max(res_img(:));
    
    bits = 2^16 - 1;
    t_new = round(rescale(res_img,0,bits));
else
    t_min = min(res_img(:))-10;
    t_new = res_img-t_min;
end
fname = 'j2k_res/cube.j2k';

imwrite(uint16(t_new), fname,'CompressionRatio',j2k_r);

acc_bits = 0;
s=dir(fname);
acc_bits = acc_bits + s.bytes*8;

read_residuals = imread(fname);
read_residuals = double(read_residuals);

if do_rescale
    read_residuals = rescale(read_residuals, t_min, t_max);
else
    read_residuals = read_residuals+t_min;
end


if decorrelation
    post_scores = read_residuals;
    post_residuals_2d = post_scores' * loadings;
    
else
    post_residuals_2d = read_residuals;
end
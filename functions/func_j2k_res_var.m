function [post_residuals_2d,acc_bits,compression_time] = func_j2k_res_var(residuals_2d, j2k_r, sub_comps)
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
tic;
[loadings,~,~,~,explained] = pca(residuals_2d);


comp = length(explained) - sub_comps;

loadings = loadings';

% disp(size(explained));

scores = loadings(1:comp,:)*residuals_2d'; 

res_img = scores;


warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir j2k_res;

t_min=min(res_img(:));
t_max=max(res_img(:));

bits = 2^16 - 1;
t_new = round(rescale(res_img,0,bits));

fname = 'j2k_res/cube.j2k';

imwrite(uint16(t_new), fname,'CompressionRatio',j2k_r);
compression_time = toc;

acc_bits = 0;
s=dir(fname);
acc_bits = acc_bits + s.bytes*8;

read_residuals = imread(fname);
read_residuals = double(read_residuals);

read_residuals = rescale(read_residuals, t_min, t_max);


post_scores = read_residuals;
post_residuals_2d = post_scores' * loadings(1:comp,:);

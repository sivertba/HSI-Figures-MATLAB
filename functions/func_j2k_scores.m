function [cube_comp,acc_bits] = func_j2k_scores(scores,loadings,dr_method,j2k_r,cube_dim)
% % 
% % Input scores: scores matrix to be compressed.
% % 
% % Input loadings: matrix used for spectral encoding
% % 
% % Input j2k_r: compression ratio set by JPEG2000 stage
% % 
% % Input cube_dim: size vector of original image cube
% % 
% % Output cube_comp: Restored cube after decompression.
% % 
% % Output acc_bits: Number of bits after compression.

scores = scores';

methods = ["None3D", "PCA3D", "None2D", "PCA2D", "OTFP","NoSpat"];

if dr_method == "PCA3D"
    prepped_scores = reshape(scores, cube_dim);
elseif dr_method == "PCA2D"
    prepped_scores = spatialize_scores(scores, cube_dim);
elseif dr_method == "OTFP"
    prepped_scores = spatialize_scores(scores, cube_dim);
else
    warning("Could not find
end

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir j2k_test;

system('rm j2k_test/*');

t_min=min(prepped_scores(:));
t_max=max(prepped_scores(:));

bits = 2^16 - 1;
t_new = round(rescale(prepped_scores,0,bits));

imwrite(uint16(t_new), 'j2k_test/png_1.png');

if j2k_r == 0
    system('opj_compress -ImgDir j2k_test/ -OutFor J2K');
else
    system(join(['opj_compress -ImgDir j2k_test/ -OutFor J2K -r ' num2str(j2k_r)]));
end

system('for i in j2k_test/*.J2K; do opj_decompress -i $i -o j2k_test/$(basename "$i" .J2K)_rev.png ; done;');

acc_bits = 0;
s=dir('j2k_test/png_1.J2K');

if dr_method == "None"
    acc_bits = acc_bits + s.bytes*8;
else
    acc_bits = acc_bits + s.bytes*8;
    acc_bits = acc_bits + numel(loadings(1:comps,:))*64;
end

fname = 'j2k_test/png_1_rev.png' ;

if dr_method == "None"
    spatialized_cube_comp = rescale(imread(fname),t_min,t_max);
    X_comp = despatialize_scores(spatialized_cube_comp, cube_dim, size(scores));
    cube_comp= reshape(X_comp,[Dsize(1) Dsize(2) Dsize(3)]);
else
    spatialized_cube_comp = rescale(imread(fname),t_min,t_max);
    scores_comp = transpose(despatialize_scores(spatialized_cube_comp, cube_dim, size(scores(1:comps,:)')));
    X_comp = scores_comp(1:comps,:)' * loadings(1:comps,:);
    cube_comp = reshape(X_comp,[Dsize(1) Dsize(2) Dsize(3)]);
end
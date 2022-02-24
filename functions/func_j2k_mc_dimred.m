function [post_cube,acc_bits] = func_j2k_mc_dimred(HSI,comps,dr_method,j2k_r, scene_name)
% % 
% % Input HSI: image cube to be compressed (spatial x spatial x spectral)
% % 
% % Input comps: How many components are kept for the dimred
% % 
% % Input dr_method: Which method for dimensionality reduction is used
% % 
% % Input j2k_r: compression ratio set by JPEG2000 stage
% % 
% % Input scene_name: Name of scene
% % 
% % Output cube_comp: Restored cube after decompression.
% % 
% % Output acc_bits: Number of bits after compression.

%%%%
M = HSI;


Dsize = size(M);
X_train= reshape(M,[Dsize(1)*Dsize(2) Dsize(3)]);

if dr_method == "None"
    pre_cube = M;

elseif dr_method == "OTFP"
    [scores, loadings, otfp_out] = otfp_prep(X_train, Dsize, scene_name);
    comps = min([otfp_out.maxdim comps]);
    pre_cube = reshape(scores(1:comps,:)',[Dsize(1) Dsize(2) comps]);

elseif dr_method == "PCA"
    [scores, loadings] = dimred_func(X_train,"PCA");
    pre_cube = reshape(scores(1:comps,:)',[Dsize(1) Dsize(2) comps]);

else
    error("Unrecognized dr_method");
end

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir j2k_test;

system('rm j2k_test/*');

t_min=min(pre_cube(:));
t_max=max(pre_cube(:));

bits = 2^16 - 1;
t_new = round(rescale(pre_cube,0,bits));

fname = 'j2k_test/cube.j2k';

imwrite(uint16(t_new), fname,'CompressionRatio',j2k_r);

acc_bits = 0;
s=dir(fname);
acc_bits = acc_bits + s.bytes*8;

read_cube = rescale(imread(fname), t_min, t_max);

if dr_method == "None"
    post_cube = read_cube;
   
elseif dr_method == "PCA"
    post_scores = reshape(read_cube,[Dsize(1)*Dsize(2) comps]);
    X_comp = post_scores * loadings(1:comps,:);
    post_cube = reshape(X_comp,[Dsize(1) Dsize(2) Dsize(3)]);

elseif dr_method == "OTFP"
    post_scores = reshape(read_cube,[Dsize(1)*Dsize(2) comps]);
    X_comp = post_scores * loadings(1:comps,:);
    post_cube = reshape(X_comp,[Dsize(1) Dsize(2) Dsize(3)]);
end

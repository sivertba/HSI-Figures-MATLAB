% WIP
function [cube_comp,acc_bits] = func_j2k_dimred_no_spat(HSI,comps,dr_method,j2k_r, scene_name)
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

M = HSI;

Dsize = size(M);
X_train= reshape(M,[Dsize(1)*Dsize(2) Dsize(3)]);

if dr_method == "OTFP"
    [scores, loadings, otfp_out] = otfp_prep(X_train, Dsize, scene_name);
    comps = min([otfp_out.maxdim comps]);
elseif dr_method == "PCA"
    [scores, loadings] = dimred_func(X_train,"PCA");
    comps = min([size(HSI,3) comps]);
else
    error("Unrecognized dr_method");
end

non_spatialized_cube = scores(1:comps,:);

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir j2k_test;

system('rm j2k_test/*');

t_min=min(non_spatialized_cube(:));
t_max=max(non_spatialized_cube(:));

bits = 2^16 - 1;
t_new = round(rescale(non_spatialized_cube,0,bits));

imwrite(uint16(t_new), 'j2k_test/png_1.png');

if j2k_r == 0
    system('opj_compress -ImgDir j2k_test/ -OutFor J2K');
else
    system(join(['opj_compress -n 1 -ImgDir j2k_test/ -OutFor J2K -r ' num2str(j2k_r)]));
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

fname = join(['j2k_test/png_1_rev.png']) ;

if dr_method == "None"
    spatialized_cube_comp = rescale(imread(fname),t_min,t_max);
    X_comp = despatialize_scores(spatialized_cube_comp, size(M), size(X_train));
elseif dr_method == "PCA"
    scores_comp = rescale(imread(fname),t_min,t_max);
    X_comp = scores_comp(1:comps,:)' * loadings(1:comps,:);
elseif dr_method == "OTFP"
    scores_comp = rescale(imread(fname), t_min, t_max);
    X_comp = scores_comp*loadings(1:comps,:);
    X_comp = X_comp ./repmat(otfp_out.total_weight,Dsize(2),1) + otfp_out.center;
end

cube_comp = reshape(X_comp,[Dsize(1) Dsize(2) Dsize(3)]);

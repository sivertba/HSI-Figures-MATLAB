function [out, ssim_vector] = MySSIM(tar,ref)
%--------------------------------------------------------------------------
% Peak Signal to Noise Ratio (PSNR)
%
% USAGE
%   out = PSNR(ref,tar)
%
% INPUT
%   ref : reference HS data (rows,cols,bands)
%   tar : target HS data (rows,cols,bands)
%
% OUTPUT
%   out.mean    : mean PSNR (scalar)
%   out.var     : variance in PSNR (scalar)
%   psnr_vector : vector of psnr value for each bandwidth
%
%--------------------------------------------------------------------------
[rows_r,cols_r,bands_r] = size(ref);
[rows_t,cols_t,bands_t] = size(tar);

assert(rows_r   == rows_t, "Wrong dimensions");
assert(cols_r   == cols_t, "Wrong dimensions");
assert(bands_r  == bands_t, "Wrong dimensions");

out = multissim3(tar,ref);

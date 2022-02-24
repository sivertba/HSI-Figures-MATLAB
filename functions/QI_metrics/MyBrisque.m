function [out, bq_vector] = MyBrisque(tar)
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
[~,~,bands_t] = size(tar);
bq_vector = zeros(bands_t,1);

for ii = 1:bands_t
    tar_2d  = tar(:,:,ii);    
    bq_vector =  brisque(tar_2d);
end

out.mean = mean(bq_vector);
out.var = var(bq_vector);

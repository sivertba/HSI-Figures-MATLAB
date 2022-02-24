function [out, snr_vector] = plainSNR(tar,ref)
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

snr_vector = zeros(bands_t,1);

for ii = 1:bands_t
    ref_2d  = ref(:,:,ii);
    tar_2d  = tar(:,:,ii);
    
    err = (norm(tar_2d(:)-ref_2d(:),2).^2)/numel(tar_2d);
    snr_vector(ii) = 10*log10(mean(ref_2d(:).^2)/err);
end

out.mean = mean(snr_vector);
out.snr_1 = out.mean;
out.snr_2 = 10 * log10( var(ref(:)) / mse(tar(:),ref(:)) );
out.var = var(snr_vector);

function out = NRMSE(tar,ref)
%--------------------------------------------------------------------------
% Root mean squared error (RMSE)
%
% USAGE
%   out = RMSE(ref,tar)
%
% INPUT
%   ref : reference HS data (rows,cols,bands)
%   tar : target HS data (rows,cols,bands)
%
% OUTPUT
%   out : NRMSE (scalar)
%
%--------------------------------------------------------------------------
% [rows,cols,bands] = size(ref);
% out = (sum(sum(sum((tar-ref).^2)))/(rows*cols*bands)).^0.5;

[rows_r,cols_r,bands_r] = size(ref);
[rows_t,cols_t,bands_t] = size(tar);

assert(rows_r   == rows_t, "Wrong dimensions");
assert(cols_r   == cols_t, "Wrong dimensions");
assert(bands_r  == bands_t, "Wrong dimensions");

nrmse_vector = zeros(bands_t,1);

for ii = 1:bands_t
    ref_2d  = ref(:,:,ii);
    tar_2d  = tar(:,:,ii);
    
    mse = (norm(tar_2d(:)-ref_2d(:),2).^2)/numel(tar_2d);
    nrmse_vector(ii) = mse /(max(max(abs(tar_2d))) - min(min(abs(tar_2d))));
end

out = mean(nrmse_vector);
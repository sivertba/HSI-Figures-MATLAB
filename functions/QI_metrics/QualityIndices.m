function Out = QualityIndices(I_HS,I_REF)
%--------------------------------------------------------------------------
% Quality Indices
%
% USAGE
%   Out = QualityIndices(I_HS,I_REF)
%
% INPUT
%   I_HS  : target HS data (rows,cols,bands)
%   I_REF : reference HS data (rows,cols,bands)
%
% OUTPUT
%   Out.cc   : CC
%   Out.sam  : SAM
%   Out.rmse : RMSE
%   Out.PSNR : mean PSNR
%--------------------------------------------------------------------------
% cc = CC(I_HS,I_REF);
% Out.cc = mean(cc);

[angle_SAM, ~] = SAM(I_HS,I_REF);
Out.sam = angle_SAM*180/pi; %degrees

% Out.nrmse = NRMSE(I_HS,I_REF);

% [psnr_val, ~] = PSNR(I_HS,I_REF);
% Out.psnr = psnr_val.mean;

[snr_val, ~] = plainSNR(I_HS,I_REF);
Out.snr = snr_val.snr_2;

% Out.ssim = MySSIM(I_HS,I_REF);

[bq_val, ~] = MyBrisque(I_HS);
Out.brisque = bq_val.mean;

return;
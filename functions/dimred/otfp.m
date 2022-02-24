%%  Computes scores and loadings with the OTFP approach
%   This is a a placeholder function as of 30/04-2020
function [scores, loadings] = otfp(X_2d)
% Input:  X_2d, A numeric matrix representing the raw spectra/radiance
%         data, with matrix dimensions (pixel, spectral response).
% Output: scores, A numeric matrix representing the scores of X_2d,
%             in the format (score vectors, score values)
% Output: loadings, A numeric matrix representing the loadings of X_2d.
%         in the format (loading vectors, loading values)

% TO DO's:
% Implement / Get OTFP code

load('initData_test.mat'); %[dataPath,'/initSatDataRun2New.mat']
loadings=initLoadingsBasis';
[scores, ~] = ProjectData(X_2d, loadings);

end

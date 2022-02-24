%%  Computes scores and loadings with Bioucas Dias' approach
function [scores, loadings] = mnf_spectral(X_2d)

% % Input:  X_2d, A numeric matrix representing the raw spectra/radiance
% %         data, with matrix dimensions (pixel, spectral response).
% % Output: scores, A numeric matrix representing the scores of X_2d,
%             in the format (score vectors, score values)
% % Output: loadings, A numeric matrix representing the loadings of X_2d.
% %         in the format (loading vectors, loading values)

% % TO DO's:
% % put comments on internal function

% compute noise vectors
[dX, ~]=estNoise(X_2d','off');

% Take the eigenvector expansion of the covariance of dX
[Un,Sn,~] = svd(dX*dX');

% Whiten the original data
Xw = X_2d*Un*inv(sqrt(Sn));

% Compute the eigenvector expansion of the covariance of wX
[~,~,Vw] = svd(Xw'*Xw);

% Define transformation matrix
loadings = Un*inv(sqrt(Sn))*Vw;
loadings = pinv(loadings);

% Compute the Maximum noise fraction scores using Green's method
scores = loadings*X_2d';
end
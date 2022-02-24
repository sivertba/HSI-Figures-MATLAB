%%  Computes scores and loadings with Green's approach
function [scores, loadings] = mnf_spatial(X_2d)

% % Input:  X_2d, A numeric matrix representing the raw spectra/radiance
% %         data, with matrix dimensions (pixel, spectral response).
% % Output: scores, A numeric matrix representing the scores of X_2d,
%             in the format (score vectors, score values)
% % Output: loadings, A numeric matrix representing the loadings of X_2d.
% %         in the format (loading vectors, loading values)

% % TO DO's:
% % Input:  homogeneus_region: A matrix representing a homogenus region of
% %         the scene to be used for estimation of the noise vectors.
% %         defaults to using the full scene for noise estimation

[N_pixels, L_responses] = size(X_2d);

% compute noise vectors
dX = zeros(N_pixels-1,L_responses);
for i=1:(N_pixels-1)
    dX(i,:) = X_2d(i,:) - X_2d(i+1,:);
end

% Take the eigenvector expansion of the covariance of dX
[Un,Sn,~] = svd(dX'*dX);

% Whiten the original data
Xw = X_2d*Un*inv(sqrt(Sn));

% Compute the eigenvector expansion of the covariance of wX
[~,~,Vw] = svd(Xw'*Xw);

% Define transformation matrix
loadings = Un*inv(sqrt(Sn))*Vw;
loadings = pinv(loadings);

% Compute the Maximum noise fraction scores using Green's method
scores = loadings*X_2d';
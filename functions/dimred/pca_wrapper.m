%%  Computes scores and loadings with of PCA in desired format
function [scores, loadings] = pca_wrapper(X_2d)
% % Input:  X_2d, A numeric matrix representing the raw spectra/radiance
% %         data, with matrix dimensions (pixel, spectral response).
% % Output: scores, A numeric matrix representing the scores of X_2d,
%             in the format (score vectors, score values)
% % Output: loadings, A numeric matrix representing the loadings of X_2d.
% %         in the format (loading vectors, loading values)

loadings_pca = pca(X_2d); 
loadings = loadings_pca(:,:)';

scores = loadings*X_2d';

end


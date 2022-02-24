%% Project data on the known model
function [scores, residuals] = ProjectData(preprocessedData, loadings)
% Input: preprocessedData, A numeric matrix representing the preprocessed spectra data.
%   NObs x Nvar
% Input: loadings, A numeric matrix representing the loading coefficients.
%   Ncomp x Nvar
% Output: scroes, A numeric matrix representing score coefficients.
%   NObs x Ncomp
% Output: residuals, A numeric matrix representing the residuals.
%   Nobs x Nvar

    % Verify inputs
    assert(nargin >= 1 && ~isempty(preprocessedData) && ...
        isnumeric(preprocessedData) && ismatrix(preprocessedData));
    assert(nargin >= 2 && ~isempty(loadings) && isnumeric(loadings) && ismatrix(loadings));
    assert(size(loadings, 2) == size(preprocessedData, 2))
    
    scores = preprocessedData * pinv(loadings);
    residuals = preprocessedData - scores * loadings;
end

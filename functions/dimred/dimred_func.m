function [scores,loadings] = dimred_func(preprocessedData, dr_method)
%--------------------------------------------------------------------------
% Compute loadings by selected bi-linear modelling method
% 
% INPUT
% preprocessedData : The centered and scaled intensity data ready for 
%                    projection.
% idleData         : The centered and scaled intensity data ready for 
%                    projection, and other OTFP data.
% 
% dr_method        : Method used for bi-linear modelling
%                      - PCA  : Plain PCA
%                      - MNFG : Green's method
%                      - MNFB : Bioucas-Diaz's method
%                      - OTFP : IDLEtechs method
%
% OUTPUT
% scores   : A numeric matrix representing the scores of preprocessedData,
%            in the format (score vectors, score values)
% 
% loadings : A numeric matrix representing the loadings of preprocessedData.
%            in the format (loading vectors, loading values)
%
%--------------------------------------------------------------------------

switch dr_method
    case 'PCA'
        [scores, loadings] = pca_wrapper(preprocessedData);
    case 'MNFG'
        [scores, loadings] = mnf_spatial(preprocessedData);
    case 'MNFB'
        [scores, loadings] = mnf_spectral(preprocessedData);
    case 'OTFP'
        [scores, loadings] = otfp(preprocessedData);    
    otherwise
        error('Unkown method for bi-linear modelling')
end

return;
end


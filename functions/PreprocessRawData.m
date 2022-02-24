%% Preprocess raw data with respect to center and weight matrices on the satellite side
function z = PreprocessRawData(inpData, centers, totalWeights, nScans, scalingMethod)
% Input: inpData, A numeric matrix representing the raw spectra data.
% Input: centers, A numeric row vector representing the center.
% Input: totalWeights: A numeric matrix representing the total weight
% Input: scalingMethod, A string representing the scaling method.
% Currently, only 'TotalWeightScaling' is supported.
% Output: z, The centered and scaled intensity data ready for projection.

    % Verify the raw data matrix input
    assert(nargin >= 1 && ~isempty(inpData) && isnumeric(inpData) && ismatrix(inpData));
    assert(nargin >= 2 && ~isempty(centers) && isnumeric(centers) && isvector(centers));
    assert(all(size(centers, 1) == 1 & size(centers, 2) == size(inpData, 2)));
    assert(nargin >= 3 && ~isempty(totalWeights) && isnumeric(totalWeights) && ismatrix(totalWeights));
    assert(size(totalWeights,2) == size(inpData,2))
    
    % verify the number of scans in the input image
    if nargin < 4
        if (size(totalWeights) == size(inpData))
            nScans = 1;
        else
            nScans=size(inpData,1)/size(totalWeights,1);
            assert(nScans==floor(nScans));
        end
    else
        assert(nargin >= 4 && ~isempty(inpData) && isnumeric(nScans));
        assert(nScans == size(inpData,1)/size(totalWeights,1))
    end
    
    % Verify the scaling method
    supportedScalingMethods = {'TotalWeightScaling'};
    if nargin < 5 || isempty(scalingMethod)
        scalingMethod = 'TotalWeightScaling';
    elseif ~ischar(scalingMethod)
        error('Invalid scaling method');
    else
        scalingMethodSupported = false;
        for methodInd = 1:length(supportedScalingMethods)
            if strcmpi(supportedScalingMethods(methodInd), scalingMethod)
                scalingMethodSupported = true;
                break;
            end
        end
        
        if ~scalingMethodSupported
            error('The scaling method is not supported');
        end
    end
        
    % Implement TotalWeightScaling
    if strcmp('TotalWeightScaling', scalingMethod)
        totalWeightsN=repmat(totalWeights,nScans,1);
        z = (inpData - centers) .* totalWeightsN;
        return;
    end
    
    error('The scaling method is supported but not implemented');
end

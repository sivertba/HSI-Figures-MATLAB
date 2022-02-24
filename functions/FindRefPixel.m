%% Find out the reference pixel in the data set
function whitePixInd = FindRefPixel(radiance)
% Input: radiance, a numeric matrix representing the radiance for all
% pixels

    % Verify inputs
    assert(nargin >= 1 && ~isempty(radiance) && isnumeric(radiance) && ismatrix(radiance));
    
    % Sort pixels by energy in a descending order
    sortedEnergy = sortrows([(1:size(radiance, 1))', sum(radiance, 2)], 2, 'descend');
    whitePixInd = sortedEnergy(1:20, 1);
    
end

%% Filter scores in a band pass fashion
function [filteredScores, filteredResiduals, filterIndex] = FilterScores(scores, residuals, scoreThresholds)
% % Input: scores, A numeric matrix representing the score matrix.
% % Input: residuals, A numeric matrix representing the residual matrix.
% % Input: seaLabels, A boolean vector representing the sea areas.
% % Input: scoreThresholds: A 3xV numeric vector representing score
% % thresholds. The first and second rows represent the upper and lower
% % thresholds for scores. The third rows represent the threshold methods.
% % Output: filteredScroes, A numeric matrix representing the scores for sea areas.
% % Output: filteredResiduals, A numeric matrix representing the residuals for sea areas.



    % Verify inputs
    assert(nargin >= 1 && ~isempty(scores) && isnumeric(scores) && ismatrix(scores));
    assert(nargin >= 2 && ~isempty(residuals) && isnumeric(residuals) && ismatrix(residuals));
    assert(size(scores, 1) == size(residuals, 1));
    assert(nargin >= 3 && ~isempty(scoreThresholds) && isnumeric(scoreThresholds));
    %assert(size(scoreThresholds, 1) == size(scores, 1));
    %assert(size(scoreThresholds, 2) == size(scores, 2));
    %assert(size(scoreThresholds, 3) == 3);
    assert(size(scoreThresholds, 1) == 3);
    assert(size(scoreThresholds, 2) == size(scores, 2));
    
    % Filter scores
    filterIndex = true(size(scores, 1), 1);
    %filterLabels(seaLabels) = true;
    for compInd = 1:size(scores, 2)
        % skip zeros
        if scoreThresholds(3, compInd) == 0
            continue
            
        % Lower-pass, method = -1
        elseif scoreThresholds(3,compInd) == -1
            filterIndex(scores(:, compInd) > scoreThresholds(1, compInd)) = false;
        
        % Higher-pass, method = 1
        elseif scoreThresholds(3,compInd) == 1
            filterIndex(scores(:, compInd) < scoreThresholds(2, compInd)) = false;
        
        % Center-pass, method = 2
        elseif scoreThresholds(3,compInd) == 2 
            filterIndex( (scores(:, compInd) < scoreThresholds(1, compInd) | ...
                scores(:, compInd) > scoreThresholds(2, compInd)) ) = false;
        
        % Out-pass, method = 3
        elseif scoreThresholds(3,compInd) == 3
            filterIndex((scores(:, compInd) > scoreThresholds(1, compInd) & ...
                scores(:, compInd) < scoreThresholds(2, compInd))) = false;
        end
    end
    
    filteredScores = scores(filterIndex, :);
    filteredResiduals = residuals(filterIndex, :);
end

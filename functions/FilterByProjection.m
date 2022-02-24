function keepIndex = FilterByProjection(X,filterSpectra,filterType)

%     wpi=FindRefPixel(X);
%     I0=mean(X(wpi,:));
%     X=X./I0;
    
    keepIndex = true(size(X,1),1);

    % Project on loading
    [pScores,pScoresResiduals] = ProjectData(X, filterSpectra);
    
    if strcmp(filterType,'SeaFilter')
        % set sea filter limit
        % lower limit, upper limit, method: -1=low pass.
        filterParam=[FindSeaFilterLimit(pScores) 0 -1]'; 
        
    else
        disp('filter type not recognized ... ')
        return
    end
    
    % Filtering on scores
    [ ~, ~,keepIndex] = FilterScores(pScores, pScoresResiduals, filterParam);

    
end
   


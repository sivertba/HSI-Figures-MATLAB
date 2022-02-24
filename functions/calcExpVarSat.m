function EV=calcExpVarSat(preprocessedData,seaScores,loadings,resIndex,sigSeaResiduals)

    mse0 = sum((preprocessedData - mean(preprocessedData) ).^2);
    recData=seaScores*loadings;
    
    if nargin ==5
        recData(resIndex,:)=recData(resIndex,:)+sigSeaResiduals;
    end
    
    mse = sum((preprocessedData - recData).^2);
    EV=100*(1-mse/mse0);
end
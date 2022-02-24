close all;
clear;
clc;
% path(pathdef);

%% path setup

hicoPath = './Data/HICO'; % image data
dataPath = './Data';
disp('paths set ...')

%% Load init data
% initial data
load([dataPath,'/initData_test.mat']) %[dataPath,'/initSatDataRun2New.mat']
load([dataPath,'/initParameters.mat'])

%% Data parameters
loadings=initLoadingsBasis';
scoreThresh=initScoreThresh;
totalWeight = initTotalWeight;

%% Load HS input
% HS input image
frameName = 'LEtestAll.mat';
load([hicoPath,'/',frameName])

% choose which scans to use from datacube (nPixels x nScans x nChannels)
firstScan = 1;
nScans = 5000;

rad = rad(:,firstScan:firstScan+nScans - 1,:);
%rgbImg = rgbImg(:,firstScan:firstScan+nScans - 1,:);

Dsize = size(rad);
X = reshape(rad,[Dsize(1)*Dsize(2) Dsize(3)]);
Nindex=1:(Dsize(1)*Dsize(2));

disp('Data and parameters loaded ... ')

%% Pre-processing

% initialize useIndex as true
useIndex=true(size(X(:,1)));

% Find I0 from input image and update weight matrix
if initParameters.reflectanceFlag
    wpi = FindRefPixel(X);
    I0 = mean(X(wpi,:));
    totalWeight = initTotalWeight./I0;
else
    I0 = ones(1,size(X,2));
end

% project on spectra representing land to identify sea pixels and filter out land pixels
if initParameters.preSeaProjectionFlag
    seaIndex = FilterByProjection(X./I0,landLoading,'SeaFilter');
    useIndex = useIndex & seaIndex;
end


if initParameters.findCenterFlag % only if using a local center, else model center is used
    % setting center as mean of the identified sea pixels
    center=mean(X(useIndex,:)./I0);
end

preprocessedData = PreprocessRawData(X, center, totalWeight, Dsize(2));

disp('Data Preprocessed ...')

%% Projection 
% Project input data onto the current loadings
[scores, residuals] = ProjectData(preprocessedData, loadings);

disp('Data Projected ...')

% Compression of scores
% !! not implemented here !!
% apply method for onboard compression of scores before calculating the
% residuals that are used in the weighted residual analysis
scores_comp = scores;

% residuals calculated using compressed scores
residuals_t_comp = preprocessedData - scores_comp * loadings;

%% Filter
% Filtering on score values based in threshold and filter method
% lower limit, upper limit, method. -1=low pass.
[filtered_scores, filtered_residuals,filterIndex] = FilterScores(scores_comp, residuals_t_comp, scoreThresh);
useIndex=useIndex & filterIndex;
filtered_scores=filtered_scores(useIndex,:);
filtered_residuals=filtered_residuals(useIndex,:);
    
disp('Scores Filtered ...')

%% Weighted residual analysis
resIndex = WeightedResidualAnalysis(filtered_residuals,200,3.5);
sig_residuals=filtered_residuals(resIndex,:);

disp('Residuals analysed...')

%% Calculate compression rate and explained variance
compressedSize = (numel(filtered_scores) + numel(useIndex) + ...
    numel(sig_residuals) + numel(resIndex)) / numel(X(useIndex,:));

compressedSizewoRes = (numel(filtered_scores) + numel(useIndex)) / numel(X(useIndex,:));

% calculated for testing, not necessary when running loop
EV=calcExpVarSat(preprocessedData(useIndex,:),filtered_scores,loadings,resIndex,sig_residuals);
EVwoRes=calcExpVarSat(preprocessedData(useIndex,:),filtered_scores,loadings);


%% Data downlinked
%seaScores,seaIndex
%sigSeaResiduals,resIndex

fprintf(' HSI input: %s\n Compressed Size = %.2f \n Compressed Size w/o residuals = %.2f \n Explained Variance =  %.4f %%\n Explained Variance w/o residuals =  %.4f %%\n Num components: %d \n Num residuals: %d (%.2f %% of pixels) \n ', frameName, compressedSize, compressedSizewoRes, EV, EVwoRes, size(loadings,1), sum(resIndex), (sum(resIndex)/(sum(useIndex)))*100);

disp('Satellite module done ...')

%% Ground

% only when reconstructing C_hat from an orthogonalized model basis
% reconstructing C_hat, estimate of concentration of known spectra
%    C_hat=filtered_scores(:,initParameters.emscIndices)*pinv(initScoreBasis);

% reconstructing image data
Nsea = find(useIndex);
Nres = Nsea(resIndex);

recData = zeros(size(useIndex,1),size(loadings,2));
recData(Nsea,:) = filtered_scores*loadings;
recData(Nres,:)=recData(Nres,:) + sig_residuals;

recData=recData./repmat(totalWeight,Dsize(2),1) + center; % Dsize(2) number of lines in the input image assumed to be known

disp('Data reconstructed from downlinked scores and residuals ...')

%% This is WIP

%% init
clc; close all; clear;

% assuming the main file is at root of the repo tree
% does not work if only running sections of code
% need to run full main.m script
path = mfilename('fullpath');
rem_str = length(mfilename());
path = path(1:end-rem_str);

all_paths = genpath(path);
addpath(all_paths);

disp('paths set ...')

meta_results_struct = struct;

%% Load init data
% initial data
load('initData_test.mat') %[dataPath,'/initSatDataRun2New.mat']
load('initParameters.mat')
load('parameters_cell.mat')

%% Load HS input
% HS input image
% frameName = 'LEtestAll.mat';
% load([hicoPath,'/',frameName])
scene = "HICO";
load('Data/data_set.mat')

data_train = data_set.hico.train;
data_test = data_set.hico.test;

% choose which scans to use from datacube (nPixels x nScans x nChannels)
firstScan = 1;
nScans = 500;

rad = data_test(:,firstScan:firstScan+nScans - 1,1:87);
%rgbImg = rgbImg(:,firstScan:firstScan+nScans - 1,:);

Dsize = size(rad);
X = reshape(rad,[Dsize(1)*Dsize(2) Dsize(3)]);
Nindex=1:(Dsize(1)*Dsize(2));

disp('Data and parameters loaded ... ')

%% Data parameters
% scoreThresh=initScoreThresh;
scoreThresh = zeros(3,20);
% totalWeight = initTotalWeight;
totalWeight = ones(500,87);
% initParameters.preSeaProjectionFlag = false;


%% Pre-processing

% initialize useIndex as true
useIndex=true(size(X(:,1)));

% Find I0 from input image and update weight matrix
if initParameters.reflectanceFlag
    wpi = FindRefPixel(X);
    I0 = mean(X(wpi,:));
    totalWeight = totalWeight./I0;
else
    I0 = ones(1,size(X,2));
end

if initParameters.findCenterFlag % only if using a local center, else model center is used
    % setting center as mean of the identified sea pixels
    center=mean(X(useIndex,:)./I0);
end

preprocessedData = PreprocessRawData(X, center, totalWeight, Dsize(2));

disp('Data Preprocessed ...')

%% Compute loadings
dr_method = 'MNFB';
[~,loadings] = dimred_func(preprocessedData,dr_method);
% [~,loadings] = dimred_func(preprocessedData,'PCA');
% [~,loadings] = dimred_func(preprocessedData,'MNFB');
% [~,loadings] = dimred_func(preprocessedData,'MNFG'); % do not use
% [~,loadings] = dimred_func(preprocessedData,'OTFP');

loadings = loadings(1:20,:);

%% Projection
% Project input data onto the current loadings
[scores, residuals] = ProjectData(preprocessedData, loadings);

disp('Data Projected ...')

res_shift = length(parameters_cell);
for param_set = 1:length(parameters_cell)
    disp(['Lossy Iteration: ' num2str(param_set)...
        '/' num2str(length(parameters_cell))]);
    
    %% Compression of scores and residuals
    compression_parameters = parameters_cell(param_set);
    compression_parameters= compression_parameters{1};
    
    [scores_comp, compressed_savings] = ...
        compress_scores(scores, compression_parameters,Dsize);
    
    % residuals calculated using compressed scores
    residuals_comp = preprocessedData - scores_comp * loadings;    
    
    %% Filter
    
    % before compression
    % Filtering on score values based in threshold and filter method
    % lower limit, upper limit, method. -1=low pass.
    [filtered_scores, filtered_residuals,filterIndex] = FilterScores(scores, residuals, scoreThresh);
    useIndex=useIndex & filterIndex;
    filtered_scores=filtered_scores(useIndex,:);
    filtered_residuals=filtered_residuals(useIndex,:);
    
    disp('Scores Filtered ...')
    
    % after compression
    % Filtering on score values based in threshold and filter method
    % lower limit, upper limit, method. -1=low pass.
    [filtered_scores_comp, filtered_residuals_comp,filterIndex_comp] = FilterScores(scores_comp, residuals_comp, scoreThresh);
    useIndex=useIndex & filterIndex;
    filtered_scores_comp=filtered_scores_comp(useIndex,:);
    filtered_residuals_comp=filtered_residuals_comp(useIndex,:);
    
    disp('Compressed Scores Filtered ...')
    
    
    %% Weighted residual analysis
    
    %before compression
    resIndex = WeightedResidualAnalysis(filtered_residuals,200,3.5);
    sig_residuals=filtered_residuals(resIndex,:);
    
    %after compression
    resIndex_comp = WeightedResidualAnalysis(filtered_residuals_comp,200,3.5);
    sig_residuals_comp=filtered_residuals_comp(resIndex_comp,:);
    
    disp('Residuals analysed...')
    
    %% Calculate compression rate and explained variance
    
    % before compression
    compressedSize = (numel(filtered_scores) + numel(useIndex) + ...
        numel(sig_residuals) + numel(resIndex)) / numel(X(useIndex,:));
    
    compressedSizewoRes = (numel(filtered_scores) + numel(useIndex)) / numel(X(useIndex,:));
    
    % calculated for testing, not necessary when running loop
    EV=calcExpVarSat(preprocessedData(useIndex,:),filtered_scores,loadings,resIndex,sig_residuals);
    EVwoRes=calcExpVarSat(preprocessedData(useIndex,:),filtered_scores,loadings);
    
    % after compression
    compressedSize_comp = ((numel(filtered_scores_comp) + numel(useIndex))*compressed_savings + ...
        numel(sig_residuals) + numel(resIndex_comp)) / numel(X(useIndex,:));
    
    compressedSizewoRes_comp = (numel(filtered_scores_comp) + numel(useIndex))*compressed_savings / numel(X(useIndex,:));
    
    % calculated for testing, not necessary when running loop
    EV_comp=calcExpVarSat(preprocessedData(useIndex,:),filtered_scores_comp,loadings,resIndex,sig_residuals);
    EVwoRes_comp=calcExpVarSat(preprocessedData(useIndex,:),filtered_scores_comp,loadings);    
    
    % Store results for plotting
    analysis_results = struct;
    analysis_results.compressedSize = compressedSize;
    analysis_results.compressedSizewoRes = compressedSizewoRes;
    analysis_results.EV = EV;
    analysis_results.EVwoRes = EVwoRes;
    analysis_results.compressedSize_comp = compressedSize_comp;
    analysis_results.compressedSizewoRes_comp = compressedSizewoRes_comp;
    analysis_results.EV_comp = EV_comp;
    analysis_results.EVwoRes_comp = EVwoRes_comp;
    
    recData = filtered_scores*loadings;
    recData=recData./repmat(totalWeight,Dsize(2),1) + center;
    recData_res = recData;
    recData_res(resIndex,:) = recData(resIndex,:) + sig_residuals;
    
    recData_comp = filtered_scores_comp*loadings;
    recData_comp=recData_comp./repmat(totalWeight,Dsize(2),1) + center;
    recData_comp_res = recData_comp;
    recData_comp_res(resIndex_comp,:) = recData_comp(resIndex_comp,:) + sig_residuals_comp;
    
    % Make Image cubes
    
    % Without residuals
    rad_ORIG = reshape(rad,[Dsize(1), Dsize(2), Dsize(3)]);
    rad_LOAD = reshape(recData,[Dsize(1), Dsize(2), Dsize(3)]);
    rad_COMP = reshape(recData_comp,[Dsize(1), Dsize(2), Dsize(3)]);

    % With residuals
    rad_LOAD_res = reshape(recData_res,[Dsize(1), Dsize(2), Dsize(3)]);
    rad_COMP_res = reshape(recData_comp_res,[Dsize(1), Dsize(2), Dsize(3)]);

    [qi_comp, qi_load] = plot_results_func(rad_ORIG,rad_LOAD, rad_COMP, ...
        analysis_results,compression_parameters,param_set, scene, dr_method);

    [qi_comp_res, qi_load_res] = plot_results_func(rad_ORIG,rad_LOAD_res, ...
        rad_COMP_res, analysis_results,compression_parameters,...
        param_set+res_shift, scene, dr_method);
    
    % Make meta struct
    meta_results_struct.("id_"+num2str(param_set)) = struct;

    meta_results_struct.("id_"+num2str(param_set)).("qi") = qi_comp;
    meta_results_struct.("id_"+num2str(param_set)).("qi_res") = qi_comp_res;
    meta_results_struct.("id_"+num2str(param_set)).("ar") = analysis_results;
    meta_results_struct.("id_"+num2str(param_set)).("cp") = compression_parameters;
      
    
end

%% Append no-compression results
compression_parameters.method = 'None';

meta_results_struct.("id_"+num2str(param_set+1)) = struct;
meta_results_struct.("id_"+num2str(param_set+1)).("qi") = qi_load;
meta_results_struct.("id_"+num2str(param_set+1)).("qi_res") = qi_load_res;
meta_results_struct.("id_"+num2str(param_set+1)).("ar") = analysis_results;
meta_results_struct.("id_"+num2str(param_set+1)).("cp") = compression_parameters;

%% Plot
plot_meta_results_func(meta_results_struct, scene, dr_method);
% get_top5_meta_results_func(meta_results_struct, scene, dr_method);


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

%% Load init data
% initial data
load('initData_test.mat') %[dataPath,'/initSatDataRun2New.mat']
load('initParameters.mat')
load('parameters_cell.mat')

allocated_bits = 10; % number of distinct values used in quantization

%% Load HS input
% HS input image
% frameName = 'LEtestAll.mat';
% load([hicoPath,'/',frameName])
load('Data/data_set.mat')
scenes = fields(data_set);
scenes = {scenes{4}};

scenario_counter = 1;

for s_num = 1:length(scenes)
    
    scene = scenes{s_num};
    data_train = data_set.(scenes{s_num}).train;
    data_test = data_set.(scenes{s_num}).test;
    
    % choose which scans to use from datacube (nPixels x nScans x nChannels)
    firstScan = 1;
    nScans = 500;
    
    bit_depth = 16;
    rad = data_test;
    start_bits = numel(rad)*bit_depth;
    
    Dsize = size(data_train);
    X_train= reshape(data_train,[Dsize(1)*Dsize(2) Dsize(3)]);
    
    Dsize = size(rad);
    X_test = reshape(rad,[Dsize(1)*Dsize(2) Dsize(3)]);
    Nindex=1:(Dsize(1)*Dsize(2));
    
    X_train = X_test; % remove later
    disp('Data and parameters loaded ... ')
    
    %% Data parameters
    % scoreThresh=initScoreThresh;
    % scoreThresh = zeros(3,20);
    totalWeight = ones(Dsize(1),Dsize(3));
    % initParameters.preSeaProjectionFlag = false;
    
    
    %% Pre-processing
    
    % initialize useIndex as true
    useIndex=true(size(X_test(:,1)));
    
    % Find I0 from input image and update weight matrix
    initParameters.reflectanceFlag = 0;
    if initParameters.reflectanceFlag
        wpi = FindRefPixel(X_test);
        I0 = mean(X_test(wpi,:));
        totalWeight = totalWeight./I0;
        
        wpi = FindRefPixel(X_train);
        totalWeight = totalWeight./I0;
    else
        I0 = ones(1,size(X_test,2));
    end
    
    if initParameters.findCenterFlag % only if using a local center, else model center is used
        % setting center as mean of the identified sea pixels
        center_test=mean(X_test(useIndex,:)./I0);
        center_train=mean(X_train(useIndex,:)./I0);
    end
    
    preprocessedData = PreprocessRawData(X_test, center_test, totalWeight, Dsize(2));
    preprocessedData_train = PreprocessRawData(X_train, center_train, totalWeight, Dsize(2));
    
    disp('Data Preprocessed ...')
    
    %% Compute loadings
    
    dr_methods = "PCA";
    dr_loadings = 50; %subspace_id(preprocessedData_train);
    
    meta_results_struct = struct;
    dr_method = "PCA";
    [~,loadings_full] = dimred_func(preprocessedData_train,dr_method);
    
    for lod_ii = 1:length(dr_loadings)
        dr_loading = dr_loadings(lod_ii);
        loadings = loadings_full(1:dr_loading,:);
        scoreThresh = zeros(3,dr_loadings(lod_ii));
        
        %% Projection
        % Project input data onto the current loadings
        [scores, ~] = ProjectData(preprocessedData, loadings);
        
        disp('Data Projected ...')
        
        res_shift = length(parameters_cell);
        score_size = size(scores);
        
        %% Compression in MATLAB
%         spatialized_scores = spatialize_scores(scores, Dsize);
%         
%         compression_parameters = parameters_cell(1);
%         compression_parameters= compression_parameters{1};
%         
%         wname = compression_parameters.wname;
%         wtreshold = compression_parameters.wtreshold;
%         
%         lvls_scores = wmaxlev(size(spatialized_scores),wname);
%         [wdec_og, lvl_book_scores] = wavedec2(spatialized_scores,lvls_scores,wname);
%         
%         keepapp = 0;
%         [~,wdec,~,~,~] = ...
%             wdencmp('gbl',wdec_og,lvl_book_scores,wname,lvls_scores,wtreshold,'h',keepapp);
%         
%         % Quantization of scores
%         q_scores_data = j2k_quantization(wdec,lvl_book_scores,allocated_bits);
%         end_bits_scores = pyHuff2bits(q_scores_data);
%         dq_scores_data = j2k_dequantization(q_scores_data,lvl_book_scores,allocated_bits);
%         scores_comp_spat =  waverec2(dq_scores_data,lvl_book_scores,compression_parameters.wname);
        
        
        %% Compression in openJpeg
        
        Ms = reshape(scores,[Dsize(1), Dsize(2), dr_loading]);
        Mx = im2uint16(Ms);
        
        system('rm j2k_test/*');
        
        for ii = 1:dr_loading
            imwrite(Mx(:,:,ii), join(['j2k_test/png_' num2str(ii) '.png']) );
        end
        
        system('opj_compress -ImgDir j2k_test/ -OutFor J2K');
        system('for i in j2k_test/*.J2K; do opj_decompress -i $i -o j2k_test/$(basename "$i" .J2K)_rev.png ; done;')
        
        acc_bits = 0;
        for ii = 1:dr_loading
            s=dir(join(['j2k_test/png_' num2str(ii) '.J2K']) );
            acc_bits = acc_bits + s.bytes*8;
        end
        
        Msad = Ms;
        for ii = 1:dr_loading
            fname = join(['j2k_test/png_' num2str(ii) '_rev.png']) ;
            Msad(:,:,ii) = im2double(imread(fname));
        end
                
        end_bits_scores2 = acc_bits;
        
%         scores_comp = despatialize_scores(scores_comp_spat, Dsize, score_size);
        scores_comp = reshape(Msad,[Dsize(1)*Dsize(2), dr_loading]);
        
        residuals = preprocessedData - scores_comp*loadings;
        
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
        [filtered_scores_comp, filtered_residuals_comp,filterIndex_comp] = FilterScores(scores_comp, residuals, scoreThresh);
        useIndex=useIndex & filterIndex;
        filtered_scores_comp=filtered_scores_comp(useIndex,:);
        filtered_residuals_comp=filtered_residuals_comp(useIndex,:);
        
        disp('Compressed Scores Filtered ...')
        
        
        %% Weighted residual analysis
        
        %before compression
        resIndex = WeightedResidualAnalysis(filtered_residuals,200,3);
        sig_residuals=filtered_residuals(resIndex,:);
        
        %after compression
        resIndex_comp = WeightedResidualAnalysis(filtered_residuals_comp,200,3);
        sig_residuals_alpha_comp=filtered_residuals_comp(resIndex_comp,:);
        
        disp('Residuals analysed...')
        
        %% Compression of residuals
        wtreshold = compression_parameters.wtreshold;
        lvls_res = wmaxlev(size(sig_residuals_alpha_comp),wname);
        [wdec_res, lvl_book_res] = wavedec2(sig_residuals_alpha_comp,lvls_res,wname);
        
        keepapp = 0;
        [~,c_res_new,~,~,~] = ...
            wdencmp('gbl',wdec_res,lvl_book_res,wname,lvls_res,wtreshold,'h',keepapp);
        
        % Quantization of residuals
        
        q_res_data = j2k_quantization(c_res_new,lvl_book_res,allocated_bits);
        end_bits_res = pyHuff2bits(q_res_data);
        dq_res_data = j2k_dequantization(q_res_data,lvl_book_res,allocated_bits);
        
        sig_residuals_comp = waverec2(dq_res_data,lvl_book_res,wname);
        
        
        %% Restoration
        
        recData = filtered_scores*loadings;
        recData=recData./repmat(totalWeight,Dsize(2),1) + center_test;
        recData_res = recData;
        recData_res(resIndex,:) = recData(resIndex,:) + sig_residuals;
        
        recData_comp = filtered_scores_comp*loadings;
        recData_comp=recData_comp./repmat(totalWeight,Dsize(2),1) + center_test;
        recData_comp_res = recData_comp;
        recData_comp_res(resIndex_comp,:) = recData_comp(resIndex_comp,:) + sig_residuals_comp;
        
        quant_err_scores = struct;
        quant_err_scores.mean = mean(c_res_new(:) - dq_res_data(:));
        quant_err_scores.var = var(c_res_new(:) - dq_res_data(:));
        
        quant_err_res = struct;
        quant_err_res.mean = mean(wdec(:) - dq_scores_data(:));
        quant_err_res.var = var(wdec(:) - dq_scores_data(:));
        
        %% Make Image cubes
        
        % Without residuals
        rad_ORIG = reshape(rad,[Dsize(1), Dsize(2), Dsize(3)]);
        rad_LOAD = reshape(recData,[Dsize(1), Dsize(2), Dsize(3)]);
        rad_COMP = reshape(recData_comp,[Dsize(1), Dsize(2), Dsize(3)]);
        
        % With residuals
        rad_LOAD_res = reshape(recData_res,[Dsize(1), Dsize(2), Dsize(3)]);
        rad_COMP_res = reshape(recData_comp_res,[Dsize(1), Dsize(2), Dsize(3)]);
        
        % Make meta struct
        meta_results_struct.("id_"+num2str(scenario_counter)) = struct;
        
        meta_results_struct.("id_"+num2str(scenario_counter)).("scene") = scene;
        meta_results_struct.("id_"+num2str(scenario_counter)).("dr_meth") = dr_method;
        meta_results_struct.("id_"+num2str(scenario_counter)).("dr_loadings") = dr_loading;
        
        meta_results_struct.("id_"+num2str(scenario_counter)).("wname") = wname;
        meta_results_struct.("id_"+num2str(scenario_counter)).("wtreshold") = wtreshold;
        
        meta_results_struct.("id_"+num2str(scenario_counter)).("start_bits") = start_bits;
        meta_results_struct.("id_"+num2str(scenario_counter)).("end_bits_res") = end_bits_res;
        meta_results_struct.("id_"+num2str(scenario_counter)).("end_bits_scores") = end_bits_scores;
        meta_results_struct.("id_"+num2str(scenario_counter)).("pixels") = numel(rad);
        meta_results_struct.("id_"+num2str(scenario_counter)).("cube_size") = size(rad);
        
        meta_results_struct.("id_"+num2str(scenario_counter)).("quant_err_scores") = quant_err_scores;
        meta_results_struct.("id_"+num2str(scenario_counter)).("quant_err_res") = quant_err_res;
        
        meta_results_struct.("id_"+num2str(scenario_counter)).("qi_load") = QualityIndices(rad,rad_LOAD);
        meta_results_struct.("id_"+num2str(scenario_counter)).("qi_comp") = QualityIndices(rad,rad_COMP);
        meta_results_struct.("id_"+num2str(scenario_counter)).("qi_load_res") = QualityIndices(rad,rad_LOAD_res);
        meta_results_struct.("id_"+num2str(scenario_counter)).("qi_comp_res") = QualityIndices(rad,rad_COMP_res);
        
        % Save meta results struct
        warning('off', 'MATLAB:MKDIR:DirectoryExists');
        mkdir results;
        
        fname_figure= "results/meta_results_struct_"+scene+"_"+dr_method+".mat";
        
        save(fname_figure,...
            'meta_results_struct')
        
        scenario_counter = scenario_counter +1;        
    end % loadings-loop
end % scene-loop

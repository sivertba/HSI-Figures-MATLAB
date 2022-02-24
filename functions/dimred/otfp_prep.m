function [scores, loadings, otfp_out] = otfp_prep(X_train, o_dims, scene_name)
load('Data/OTFPmodels.mat');
otfp_out = OTFPmodels.(join([lower(scene_name) , '_Test'],''));

totalWeight = ones(o_dims(1),o_dims(3));

% initialize useIndex as true
useIndex=true(size(X_train(:,1)));

% Find I0 from input image and update weight matrix
I0 = otfp_out.I0;
totalWeight = totalWeight./I0;
center = otfp_out.center;


preprocessedData = PreprocessRawData(X_train, center, totalWeight, o_dims(2));

loadings = pinv(otfp_out.loadings)';
scores = loadings*preprocessedData';

otfp_out.maxdim = size(loadings,1);
otfp_out.total_weight = totalWeight;

end

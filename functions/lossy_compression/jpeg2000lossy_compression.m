%%  jpeg2000lossy compressed scores and percentage of thresholded coefficients 
%   i.e. coefficients that are not equal to 0
function [compressed_scores, thresholded_coefficients, transform] = jpeg2000lossy_compression(original_scores, parameters)
% % Input:  original_scores, A numeric matrix representing the scores.
% % Input:  parameters, A struct structure with parameters for fft
% %             - method, string saying fft.
% %             - CR_goal, scalar value of target compression ratio
% % 
% % Output: compressed_scores, A numeric matrix representing the scores of 
% %         after compression with fft.
% % Output: thresholded_coefficients, A scalar representing the Compression 
% %         score, the percentage of thresholded coefficients that are not 0.
% % Output: transform, A struct containing the variables needed for the 
% %         transformation between fourier space and orgiinal


assert(parameters.method == "jpeg2000lossy", "Wrong set of parameters");

[sz1, sz2] = size(original_scores);
data_matrix_vec = reshape(original_scores,[1 sz1*sz2]);
 
% assert(reshape(data_matrix_vec, [sz1,sz2]) == original_scores, ...
%     "something wrong with reshape in function");

sz_ratio = (64/16);
data_img_uint16 = typecast(data_matrix_vec,'uint16');
data_img_uint16 = reshape(data_img_uint16, [sz1,sz2*sz_ratio]);

temp_fpath = 'endcoded.jp2';

imwrite(data_img_uint16,temp_fpath,...
    'Mode','lossy','CompressionRatio',15);

% Number of bytes originally
bytes_before = whos('original_scores');
bytes_before = bytes_before.bytes;

bytes_after = dir(temp_fpath);
bytes_after = bytes_after.bytes;

thresholded_coefficients = bytes_after/bytes_before;

data_img = imread(temp_fpath);

data_vec = reshape(data_img, [1 sz1*sz2*sz_ratio]);
data_vec_double = typecast(data_vec,'double');
compressed_scores = reshape(data_vec_double, [sz1,sz2]);

delete(temp_fpath);

transform = 0;
end
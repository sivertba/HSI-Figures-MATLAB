%% Compress scores with method specified in parameters.method
function [compressed_scores, thresholded_coefficients, transform] = compress_scores(original_scores, parameters, cube_size)
% % Input:  original_scores, A numeric matrix representing the scores.
% % Input:  parameters, A struct structure with parameters for a
% %         compression method (fft, wavelet, ...)
% % Input:  cube_size, A 1x3 vector representing the dimensions of the 
% %         image cube originally
% % 
% % Output: compressed_scores, A numeric matrix representing the scores of 
% %         after compression with wavelets.
% % Output: thresholded_coefficients, A scalar representing the Compression 
% %         score, the percentage of thresholded coefficients that are not 0.
% % Output: transform, A struct containing the variables needed for the 
% %         transformation between compression method space and orgiinal

method = parameters.method;

disp(['Using method ' method ' for compression of scores'])
disp(['With parameters:']);
disp(parameters);

spatialized_scores = spatialize_scores(original_scores, cube_size);

switch method
    case "fft"
        [cs, tc, tf] = fft_compression(spatialized_scores, parameters);
    case "wavelet"
        [cs, tc, tf] = wavelet_compression(spatialized_scores, parameters);
    case "jpeg2000lossy" % does not work
        [cs, tc, tf] = jpeg2000lossy_compression(spatialized_scores, parameters);    
    otherwise
        error(["Method" method " not found in supported methods"]);
end

score_size = size(original_scores);

compressed_scores = despatialize_scores(cs, cube_size, score_size);
thresholded_coefficients = tc;
transform = tf;

end
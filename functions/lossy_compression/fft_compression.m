%%  FFT compressed scores and percentage of thresholded coefficients 
%   i.e. coefficients that are not equal to 0
function [compressed_scores, thresholded_coefficients, transform] = fft_compression(original_scores, parameters)
% % Input:  original_scores, A numeric matrix representing the scores.
% % Input:  parameters, A struct structure with parameters for fft
% %             - method, string saying fft.
% %             - tresh, scalar value of dimensionless treshold.
% % 
% % Output: compressed_scores, A numeric matrix representing the scores of 
% %         after compression with fft.
% % Output: thresholded_coefficients, A scalar representing the Compression 
% %         score, the percentage of thresholded coefficients that are not 0.
% % Output: transform, A struct containing the variables needed for the 
% %         transformation between fourier space and orgiinal


assert(parameters.method == "fft", "Wrong set of parameters");

fft_tresh_dimless = parameters.tresh;
s_size = size(original_scores);

fft_scores = fft2(original_scores); 

fft_tresh = fft_tresh_dimless * max(max(abs(fft_scores)));
fft_ind = abs(fft_scores)>fft_tresh;
fft_count = s_size(1) * s_size(2) - sum(sum(fft_ind));
fftlow_scores = fft_scores.*fft_ind;

thresholded_coefficients = 1-fft_count/(s_size(1)*s_size(2));
compressed_scores = ifft2(fftlow_scores);

transform = struct;
transform.fourier_img = fftlow_scores;

end
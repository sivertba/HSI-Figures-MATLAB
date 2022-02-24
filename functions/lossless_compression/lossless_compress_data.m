%% Compress scores with method specified in parameters.method
function bytes_used = lossless_compress_data(transform_struct, method)
% % Input:  original_data, A struct representing the transform
% % Input:  parameters, A struct structure with parameters for a
% %         lossy compression method used prior
% % Input:  method, A struct structure with parameters for a
% %         lossless compression method e.g. source encoding
% % 
% % Output: bytes_used, Number of bytes used in new representation
% % Output: transform, A struct containing the variables needed for the 
% %         transformation between compression method space and orgiinal

switch nargin
    case 1
        disp('No method selected, will use all and do votation');
    case 2
        disp(method);
    otherwise
        error("Wrong number of paramters");
end

switch method
    case "rle"
        bytes_used = rle_compression(transform_struct);
    case "jpeg2000"
        bytes_used = jpeg2000_lossless_compression(transform_struct);
    otherwise
        error(["Method" method " not found in supported methods"]);
end

compressed_scores = cs;
thresholded_coefficients = tc;
transform = tf;

imwrite
end
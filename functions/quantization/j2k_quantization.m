%%  Quantization of wavelet bitplanes.
function q_data = j2k_quantization(data,lvl_book,allocated_bits)
% % Input:  data, a vector cointaining wavlet values to be quantized
% % Input:  lvl_book, Bookkeeping matrix. The matrix S contains 
% %                   the dimensions of the wavelet coefficients by level 
% %                   and is used to parse the wavelet values
% % Input:  allocated_bits, max number of distinct values from quantization
% %                         denoted as log2(allocated_bits)
% % 
% % Output: q_data, data quantized to the desired granulation


R_b = 0; % this values does not seem to do anything ...


q_data = data;
cursor = 1;
eta_0 = allocated_bits;
mu_b = allocated_bits;


for bp = 2:size(lvl_book,1)-1
    if bp == 2
        l_bp = 4*lvl_book(bp,1)* lvl_book(bp,2);
    else
        l_bp = 3*lvl_book(bp,1)* lvl_book(bp,2);        
    end
        
    n_b = bp-2;
    eta_b = eta_0 + n_b - (size(lvl_book,1) - bp);
    
    delta_b = 2^(R_b-eta_b)*(1+mu_b/2^11);

    for bp_ii = cursor:cursor+l_bp-1
        q_data(bp_ii) = ...
            sign(data(bp_ii))*floor(abs(data(bp_ii))/delta_b);
    end
     
    cursor = cursor +l_bp;
end

end


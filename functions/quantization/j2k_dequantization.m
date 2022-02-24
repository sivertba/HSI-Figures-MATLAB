%%  Quantization of wavelet bitplanes.
function dq_data = j2k_dequantization(q_data,lvl_book,allocated_bits)
% % Input:  data, a vector cointaining quantized wavlet values
% % Input:  lvl_book, Bookkeeping matrix. The matrix S contains 
% %                   the dimensions of the wavelet coefficients by level 
% %                   and is used to parse the wavelet values
% % Input:  allocated_bits, max number of distinct values from quantization
% %                         denoted as log2(allocated_bits)
% % 
% % Output: dq_data, de-quantized  data

eta_0 = allocated_bits;
mu_b = allocated_bits;
R_b = 0;

dq_data = q_data;
cursor = 1;
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
        if q_data(bp_ii) > 0
            dq_data(bp_ii) = ...
                q_data(bp_ii)*delta_b;
        
        elseif q_data(bp_ii) < 0
            dq_data(bp_ii) = ...
                q_data(bp_ii)*delta_b;
        
        elseif q_data(bp_ii) == 0
            dq_data(bp_ii) = 0;
            
        end
    end
     
    cursor = cursor +l_bp;
end


end


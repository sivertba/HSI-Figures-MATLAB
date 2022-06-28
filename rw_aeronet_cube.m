clc; close all; clear;
%% Load

fname       = "Data/frohavet_2022-04-19.bip";
DimsHWN     = [956,684,120];
precision   = "uint16";
offset      = 0;
interleave  = "bip";
byteorder   = 'ieee-le';

cube = multibandread(fname,DimsHWN,precision,offset,interleave,byteorder);

%% Display

spectral_coeffs = [ -5.719788129534360902e-09 
                    1.324037080791479811e-05 
                    3.751455956374321055e-01
                    2.264762366937773663e+02];
x_start = 428;
x_stop = 1508;
image_width = 120;
x = linspace(x_start, x_stop, image_width);
h1_wl = spectral_coeffs(4) + ...
        spectral_coeffs(3).*x + ...
        spectral_coeffs(2).*x.*x + ...
        spectral_coeffs(1).*x.*x.*x;

[~, Red] = min(abs(h1_wl-665));
[~, Green] = min(abs(h1_wl-540));
[~, Blue] = min(abs(h1_wl-470));


%% Resize

h1_data_dim = [956,684,120];
h1_data = zeros(h1_data_dim);

for i = 1:size(cube,3)
    h1_data(:,:,i) = imresize(cube(:,:,i), h1_data_dim(1:2));
end



save("h1_data.mat","h1_data","h1_wl","h1_data_dim");



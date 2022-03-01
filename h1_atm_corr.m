clc; close all; clear;
%% Load
load('h1_data.mat')

[~, Red]    = min(abs(h1_wl-665));
[~, Green]  = min(abs(h1_wl-540));
[~, Blue]   = min(abs(h1_wl-470));

dir = "figs/atm_corr/";

h1_sets = struct;
% iarr - Does not work
h1_sets.iarr = iarr(h1_data);

% log Residuals - Does not work
h1_sets.logR = logResiduals(h1_data);

% There is no signal prior to 400 nm and this makes the ff throw an error
h1_darkpixel = subtractDarkPixel(h1_data(:,:,4:end));
h1_sets.ff = flatField(h1_darkpixel, [85 170 5 10]);

h1_fields = fields(h1_sets);
for ii = 1:length(h1_fields)
    h1 = figure;

    current_field = h1_fields{ii};
    current_set = h1_sets.(current_field);
    subplot(1,3,1);
    create_RGB(h1_data,Red,Green,Blue);
    title("Original");

    subplot(1,3,2);
    create_RGB(current_set,Red,Green,Blue);
    title(current_field);

    subplot(1,3,3);
    if current_field == "ff"
       normalized_difference = h1_data(:,:,4:end)./max(h1_data(:,:,4:end)) - current_set./max(current_set);
    else
       normalized_difference = h1_data./max(h1_data) - current_set./max(current_set);
    end
    normalized_difference = max(normalized_difference,[],3);
    imshow(normalized_difference);
    title("Normalized Max Difference")
    
    set(gcf, 'Position', get(0, 'Screensize'));
    exportgraphics(h1,sprintf("%splots_%s.png",dir,current_field));

end

%% Mask

[~, g1]     = min(abs(h1_wl-555));
[~, b1]     = min(abs(h1_wl-443));
[~, b2]     = min(abs(h1_wl-490));
[~, b3]     = min(abs(h1_wl-510));



green_nm = h1_wl(g1);
blue1_nm = h1_wl(b1);
blue2_nm = h1_wl(b2);
blue3_nm = h1_wl(b3);

ax = [0.3272, -2.9940, 2.7218, -1.2259, -0.5683];

h1_masked = h1_data;
h1_chla = zeros(size(h1_data,1), size(h1_data,2));

h1_ff = h1_sets.ff;

mask_coeff = 5;
for xx = 1:size(h1_data,1)
    for yy = 1:size(h1_data,2)

        if h1_data(xx,yy,108) >= mask_coeff
            h1_masked(xx,yy,:) = zeros(1,size(h1_data,3));
            h1_chla(xx,yy) = 0;
        else
            chlor_a = ax(1);

            Rrs_green = h1_ff(xx,yy,g1);
            Rrs_blue = max([...
                h1_ff(xx,yy,b1),...
                h1_ff(xx,yy,b2),...
                h1_ff(xx,yy,b3)]);

            for i = 2:length(ax)
                chlor_a = chlor_a + ax(i)*(log10(Rrs_blue/Rrs_green))^i;
            end
            
            chlor_a = 10^(chlor_a);
            h1_chla(xx,yy) = chlor_a;
%             disp(chlor_a);

        end
    end
end

final_fig = figure;
subplot(1,4,1)
create_RGB(h1_data,Red,Green,Blue);
title("Orignial");

subplot(1,4,2)
create_RGB(h1_ff,Red,Green,Blue);
title("FF Atm.Corr");

subplot(1,4,3)
create_RGB(h1_masked,Red,Green,Blue);
title("Land Masked");

subplot(1,4,4)
imshow(h1_chla);
title("Chl-a estimate");

set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(final_fig,sprintf("%splots_chla.png",dir));
%%
close all;

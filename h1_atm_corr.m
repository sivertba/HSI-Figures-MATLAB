clc; close all; clear;
%% Load
load('h1_data.mat')

[~, Red]    = min(abs(h1_wl-665));
[~, Green]  = min(abs(h1_wl-540));
[~, Blue]   = min(abs(h1_wl-470));

dir = "figs/atm_corr/";

h1_sets = struct;
h1_sets.original = h1_data;

% iarr - Does not work
h1_sets.iarr = iarr(h1_data);

% log Residuals - Does not work
h1_sets.logR = logResiduals(h1_data);

% There is no signal prior to 400 nm and this makes the ff throw an error
h1_darkpixel = subtractDarkPixel(h1_data(:,:,4:end));
h1_sets.ff = flatField(h1_darkpixel, [85 170 5 10]);

h1_fields = fields(h1_sets);
for ii = 1:length(h1_fields)
    current_field = h1_fields{ii};
    current_set = h1_sets.(current_field);

    h3 = figure;
    create_RGB(current_set,Red,Green,Blue);
    title(current_field);

    set(gcf, 'Position', get(0, 'Screensize'));
    exportgraphics(h3,sprintf("%splots_method_%s.png",dir,current_field));
    
    h4 = figure;
    if current_field == "ff"
       normalized_difference = h1_data(:,:,4:end)./max(h1_data(:,:,4:end)) - current_set./max(current_set);
    else
       normalized_difference = h1_data./max(h1_data) - current_set./max(current_set);
    end
    normalized_difference = max(normalized_difference,[],3);
    imagesc(normalized_difference); axis off; axis image; colorbar; colormap("turbo");
    title(sprintf("%s Max Difference", current_field));
    
    set(gcf, 'Position', get(0, 'Screensize'));
    exportgraphics(h4,sprintf("%splots_max_diff_%s.png",dir,current_field));
%     break;
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

f1 = figure;
create_RGB(h1_masked,Red,Green,Blue);
title("Land Masked");
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(f1,sprintf("%splots_masked.png",dir));

f2 = figure;
% imshow(h1_chla);
imagesc(h1_chla); axis off; axis image; colorbar; colormap("turbo");
title("Chl-a estimate");
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(f2,sprintf("%splots_chla.png",dir));

%% Plot Water Pixels in sets

Water_region_yy = 100:150; 
Water_region_xx = 600:640; 
Water_data = h1_data(Water_region_xx,Water_region_yy,:);
for ii = 1:length(h1_fields)
    current_field = h1_fields{ii};
    current_set = h1_sets.(current_field);

    Water_data = current_set(Water_region_xx,Water_region_yy,:);
    Water_data = reshape(Water_data,...
        [size(Water_data,1)*size(Water_data,2) size(Water_data,3)]);

    hx = figure;

    if current_field == "ff"
        plot(h1_wl(4:end), Water_data');
    else
        plot(h1_wl, Water_data');
    end

    xlabel("Wavelength (nm)");
    ylabel("Reflectance or Radiance");
    title(current_field + " Water Spectra")
    
    exportgraphics(hx,sprintf("%s%s_water_spectra.png",dir,current_field));


end


%%
close all;

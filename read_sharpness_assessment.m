%% spatial analysis
clc; clear; close all;

%% Load data
load("h1_data.mat")
xlsfname = 'mjosa_2022-06-23_stats.xls';
[~,sheet_name] = xlsfinfo(xlsfname);
for k=1:numel(sheet_name)
  data{k} = readtable(xlsfname,"Sheet",sheet_name{k});
end

%% Collect
num_edges = zeros(size(h1_wl));
fwhm_means = zeros(size(h1_wl));
fwhm_stds = zeros(size(h1_wl));

for ii = 1:size(data,2)
    current_band = data{ii};
    fwhm_edges = current_band.Var2(33);
    fwhm_mean = current_band.Var3(33);
    fwhm_std = current_band.Var3(34);

    num_edges(ii) = fwhm_edges;
    fwhm_means(ii) = fwhm_mean;
    fwhm_stds(ii) = fwhm_std;
end

%% Plot
hold on;

title("Sharpness Assessment for HYPSO-1 Payload");

yyaxis left;
bar(h1_wl,num_edges);
xlabel("Wavelength (nm)");
ylabel("Eligible Edges");

yyaxis right;
errorbar(h1_wl,fwhm_means,fwhm_stds, "LineWidth",1,"Marker",".", "MarkerSize",20);
ylabel("Estimated Pixel FWHM");

xlim([470 750]);

exportgraphics(gcf,"h1_sharpness.pdf",Resolution=600);

hold off;
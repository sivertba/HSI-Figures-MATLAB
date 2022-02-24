% clear all;
close all;
clc;

%% Load full image
A=imread('sample_hico_bw','png');

%% Make image black and white
B=rgb2gray(A);
[nx,ny]=size(B);

%% Compute the wavelets of our image using dwt2
n = 2;
w = 'sym2';
[C,S] = wavedec2(B,n,w);

%% Wavelet decomposition (2 level)
% LEVEL 1
A1 = appcoef2(C,S,w,1); % Approximation
[H1, V1, D1] = detcoef2('a',C,S,1); % Details
A1 = wcodemat(A1,128);
H1 = wcodemat(H1,128);
V1 = wcodemat(V1,128);
D1 = wcodemat(D1,128);
% LEVEL 2
A2 = appcoef2(C,S,w,2); % Approximation
[H2, V2, D2] = detcoef2('a',C,S,2); % Details
A2 = wcodemat(A2,128);
H2 = wcodemat(H2,128);
V2 = wcodemat(V2,128);
D2 = wcodemat(D2,128);
dec2 = [A2 H2; V2 D2];
dec1 = [imresize(dec2,size(H1)) H1 ; V1 D1];

fig1 = figure;
imshow(B);

fig2 = figure;
imshow(imresize(dec2/max(dec2(:)),size(B)));

fig3 = figure;
imshow(imresize(dec1/max(dec1(:)),size(B)));

exportgraphics(fig1,'figs/wavelet_demo_1.pdf')
exportgraphics(fig2,'figs/wavelet_demo_2.pdf')
exportgraphics(fig3,'figs/wavelet_demo_3.pdf')

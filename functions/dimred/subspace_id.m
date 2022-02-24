%% Hyperspectral signal subspace estimation
function [kf, Ek]=subspace_id(X_2d)
% % Input:  X_2d, A numeric matrix representing the raw spectra/radiance
% %         data, with matrix dimensions (pixel, spectral response).
% %        
% % Output: kf, signal subspace dimension
% % Output: Ek, matrix which columns are the eigenvectors that span 
% %             the signal subspace

assert(numel(X_2d) ~= 0,'the data set is empty');

y = X_2d';

[L, N] = size(y);

[w, Rn] = estNoise(y,'off');

[~, ~] = size(w);

x = y - w;

Ry = y*y'/N;   % sample correlation matrix 
Rx = x*x'/N;   % signal correlation matrix estimates 

[E,~]=svd(Rx); % eigen values of Rx in decreasing order, equation (15)

Rn=Rn+sum(diag(Rx))/L/10^5*eye(L);

Py = diag(E'*Ry*E); %equation (23)
Pn = diag(E'*Rn*E); %equation (24)
cost_F = -Py + 2 * Pn; %equation (22)
kf = sum(cost_F<0);
[~,ind_asc] = sort( cost_F ,'ascend');
% Ek = E(:,ind_asc(1:kf));
Ek = E(:,ind_asc(1:end));

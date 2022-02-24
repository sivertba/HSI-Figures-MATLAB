function [rxScore] = rx_anomaly(X3D,confCoefficient = 0.98;)

Osize = size(Data);
X2d = reshape(Data,[Osize(1)*Osize(2) Osize(3)]);

Cov_x = cov(X2d);

mu_x = mean(X2d,1);

RX_mat = zeros(size(X2d,1),1);

for ii = 1:length(X2d)
    r = X2d(ii,:);
    
    Rx_val = (r-mu_x) * inv(Cov_x)* (r-mu_x)';
    
    RX_mat(ii) = Rx_val;
end

%%
rxScore = reshape(RX_mat,[Osize(1) Osize(2)]);

% count = imhist(rxScore);
% pdf = count/prod(size(rxScore,[1 2]));
% cdf = cumsum(pdf(:));
% 
% confCoefficient = 0.98;
% rxThreshold = find(cdf > confCoefficient,1);
% 
% bw = rxScore > rxThreshold;

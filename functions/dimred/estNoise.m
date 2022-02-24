function [varargout]=estNoise(varargin)
narginchk(1, 3)
if nargout > 2, error('too many output parameters'); end
y = varargin{1};
if ~isnumeric(y), error('the data set must an L x N matrix'); end
noise_type = 'additive'; % default value
verbose = 1; verb ='on'; % default value
for i=2:nargin
    switch lower(varargin{i})
        case {'additive'}, noise_type = 'additive';
        case {'poisson'}, noise_type = 'poisson';
        case {'on'}, verbose = 1; verb = 'on';
        case {'off'}, verbose = 0; verb = 'off';
        otherwise, error('parameter [%d] is unknown',i);
    end
end

[L, N] = size(y);
if L<2, error('Too few bands to estimate the noise.'); end

if verbose, fprintf(1,'Noise estimates:\n'); end

if strcmp(noise_type,'poisson')
    sqy = sqrt(y.*(y>0));          % prevent negative values
    [u , ~] = estAdditiveNoise(sqy,verb); % noise estimates
    x = (sqy - u).^2;            % signal estimates
    w = sqrt(x).*u*2;
    Rw = w*w'/N;
else % additive noise
    [w, Rw] = estAdditiveNoise(y,verb); % noise estimates
end

varargout(1) = {w};
if nargout == 2, varargout(2) = {Rw}; end
return
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,Rw]=estAdditiveNoise(r,verbose)

small = 1e-6;
verbose = ~strcmpi(verbose,'off');
[L, N] = size(r);
% the noise estimation algorithm
w=zeros(L,N);
if verbose
    fprintf(1,'computing the sample correlation matrix and its inverse\n');
end
RR=r*r';     % equation (11)
RRi=inv(RR+small*eye(L)); % equation (11)
if verbose, fprintf(1,'computing band    ');end
for i=1:L
    if verbose, fprintf(1,'\b\b\b%3d',i);end
    % equation (14)
    XX = RRi - (RRi(:,i)*RRi(i,:))/RRi(i,i);
    RRa = RR(:,i); RRa(i)=0; % this remove the effects of XX(:,i)
    % equation (9)
    beta = XX * RRa; beta(i)=0; % this remove the effects of XX(i,:)
    % equation (10)
    w(i,:) = r(i,:) - beta'*r; % note that beta(i)=0 => beta(i)*r(i,:)=0
end
if verbose, fprintf(1,'\ncomputing noise correlation matrix\n');end
Rw=diag(diag(w*w'/N));
return
end
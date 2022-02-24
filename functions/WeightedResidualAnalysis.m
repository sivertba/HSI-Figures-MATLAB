%% Weighted residual analysis
function identifiedRows = WeightedResidualAnalysis(residuals,nMaxIter,nSigma)
% Input: 
%   residuals, A numeric matrix representing the residuals.
%   nMaxIter, max number of iterations of the FitGaussDist function
%   nSigma, the limit of number of standard deviations for finding significant
%   residuals. Default is set to 3.
% Output: identifiedRows, A logical vector representing the identified unmodeled structure.

    % Verify the residual matrix input
    assert(nargin >= 1 && ~isempty(residuals) && ...
        isnumeric(residuals) && ismatrix(residuals) && ...
        size(residuals, 1) >= 2);
    
    if nargin < 2
        nMaxIter=100;
    else
        assert(isnumeric(nMaxIter));
    end
    if nargin < 3
        nSigma = 3;
    else
        assert(isnumeric(nSigma));
    end


    nObs = size(residuals, 1);
    identifiedRows = false(nObs, 1);
    for bandInd = 1:size(residuals, 2)
        [counts, x] = hist(residuals(:, bandInd), 1000);
        y = counts ./ sum(counts);
        [~, estMu, estSigma] = FitGaussDist(y, x, nMaxIter);
        resLB = estMu - estSigma * nSigma;
        resUB = estMu + estSigma * nSigma;
        identifiedRows = identifiedRows | ...
            residuals(:, bandInd) < resLB | ...
            residuals(:, bandInd) > resUB;
    end
end

%% Fit a Normal Distribution to the data
function [area, mu, sigma] = FitGaussDist(yk, xk, nMaxIters)
% Input: yk, input Y sample values (1xM)
% Input: xk, input X sample values (1xM)
% Output: A, fitted area
% Output: mu, fitted mean
% Output: sigma, fitted standard deviation
    
    if nargin<3
        nMaxIters = 100;
    end
    
    Chi2RelTol = 1e-3;

    % Check input arguments
    assert(nargin >= 1 && ~isempty(yk) && isnumeric(yk) && isvector(yk));
    assert(nargin >= 2 && ~isempty(xk) && isnumeric(xk) && isvector(xk));
    assert(length(yk) == length(xk));

	% Guess mean
	idx = FindMaxima(yk);
	if isempty(idx)
		% No peak (probably a peak off the edge). Poor chance of convergence. The
		% best we can do is to start MU at highest the location of the maximum.
		[~, idx2] = max(yk);
		mu = xk(idx2);
	else
		% We have peak(s). Set to highest peak.
		[~, idx2] = max(yk(idx));
		mu = xk(idx(idx2));
    end

    % Guess sigma to be about halfway down from peak.
	mk = yk > (max(yk) / 2 + max(yk) / length(yk));
	tophalf = yk(mk); % Top half values
	st = find(yk == tophalf(1)); % Only use first value at top half
	sigma = abs(mu - xk(st(1))); % now SD=mu-1st value over FWHM
	
    % Guess for sigma might be zero !!
	if abs(sigma) < 1e-8
		% Stab in the dark, could be improved
		sigma = xk(end) / 2;
    end
    
    area = Traprule(yk, mean(diff(xk)));

    % Form A or rather Ai
    Ai = [area, mu, sigma];
    Ai = Ai(:);

    % Iterative optimization
    counter = 0;
    loopInd = 1;
    lambda  = 1000;
    fold = [0.01, 0.01, 0.01];
    fold = fold(:);
    
    while loopInd == 1
        counter = counter + 1;
        Ah = Ai; % Store old values
        area = Ai(1); % Area factor
        mu = Ai(2); % Mean
        sigma = Ai(3); % Standard deviation
        fx = GaussFunc(xk, area, mu, sigma); % abbreviation ... xk = independent variable
        rk = yk - fx;
        
        % First derivatives
        dfdA = fx ./ area;
        dfdmu = fx .* (xk - mu) ./ sigma^2;
        dfdsigma = fx .* ((xk - mu).^2 / sigma^3 - 1 / sigma);
        
        % Form Jacobian
        Jacob(1, 1) = sum((1 + lambda) * dfdA.^2);
        Jacob(1, 2) = sum(dfdA .* dfdmu);
        Jacob(1, 3) = sum(dfdA .* dfdsigma);
        Jacob(2, 1) = Jacob(1, 2);
        Jacob(2, 2) = sum((1 + lambda) * dfdmu.^2);
        Jacob(2, 3) = sum(dfdmu .* dfdsigma);
        Jacob(3, 1) = Jacob(1, 3);
        Jacob(3, 2) = Jacob(2, 3);
        Jacob(3, 3) = sum((1 + lambda) * dfdsigma.^2);
        
        % Form F
        F(1) = sum(rk .* dfdA);
        F(2) = sum(rk .* dfdmu);
        F(3) = sum(rk .* dfdsigma);
        F = F(:);
        Ai = Ah + Jacob \ F;
        rkold   = yk - GaussFunc(xk, Ah(1), Ah(2), Ah(3));  % for testing the end of the loop
        chi2old = sumsq(rkold);
        rk = yk - GaussFunc(xk, Ai(1), Ai(2), Ai(3));
        chi2 = sumsq(rk);
        
        % Calculate relative change in Chi^2 from previous iteration.
        pc_change = (sumsq(F) - sumsq(fold)) / sumsq(fold);
        fold = F;

        if chi2 > chi2old
            % Residuals have INCREASED. Make refinements coarser.
            lambda = lambda*10;
        else
            if (chi2-chi2old < 1) && (abs(pc_change) < Chi2RelTol)
                loopInd = 0;
            end
            
            % Residuals have decreased, make finer refinements.
            lambda = lambda/10;
        end
        
        if counter > nMaxIters
            loopInd = 0;
        end
        
        if isnan(chi2)
            area = NaN;
            mu = NaN;
            sigma = NaN;
            return;
        end
    end
end

%% Abbreviation to evaluate Gaussian functions.
function output = GaussFunc(x, A, mu, sigma)
    output = A * exp(-0.5 * ((x - mu) ./ sigma).^2) / sqrt(2 * pi * sigma^2);
end

%% Calculate the sum of the square of each element in the matrix
function a = sumsq(b, dim)
    if nargin == 1
        a = sum(b .* b);
        return;
    end
    
    if ndims(b) < dim
        error('DIM exceeds number of dimensions of B.');
    end
    
    a = sum(b .* b, dim);
end

%% Finds indices of y where y is locally maximum.
function idx = FindMaxima(y)
    if length(y) < 3
        idx = [];
        return;
    end
    
    grad = diff(y);
    idx = find(grad(2:end) < 0 & grad(1:end-1) > 0) + 1;
end

%% Approximate the area under a curve using the trapezium rule
function output = Traprule(y, sep)
% Input: y, y values at each x sample
% Input: sep, x separation between samples
% Output: output, calculate area
    if nargin == 2
        output = sep*(sum(y(:)) - 0.5 * (y(1) + y(end)));
    else
        output = sum(y(:)) - 0.5 * (y(1) + y(end));
    end
end

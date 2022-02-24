function out = lwt2_levels(X,wname, tresh, levels)

maxlvl = floor(log2(min(size(X))));

switch nargin
    case 1
        wname = 'haar';
        tresh = 0;
        levels = maxlvl;
    case 2
        tresh = 0;
        levels = maxlvl;
    case 3
        levels = maxlvl;
    case 4
        % do nothing
    otherwise
        warning("Wrong number of arguments");
end

assert( maxlvl >= levels,"Too many levels");

% inty
els = {'p',[-0.125 0.125],0};
lsnewInt = addlift(liftwave(wname),els);

out = X;

for ii = 1:levels
    
    t = out(1:2^(1+maxlvl-ii), 1:2^(1+maxlvl-ii));
    
    [cAint,cHint,cVint,cDint] = lwt2(t,lsnewInt);
    
    out(1:2^(1+maxlvl-ii), 1:2^(1+maxlvl-ii)) = ...
        [cAint, cHint; cVint, cDint];
    
    
end

if tresh ~= 0
    
    temp = out;
    
    for xx = 1:size(temp,1)
        for yy = 1:size(temp,2)
            tar = temp(xx,yy);
            if tar <= tresh && tar >= 0
                tar = 0;
            elseif tar >= -tresh && tar <= 0
                tar = 0;
            end
            
            temp(xx,yy) = tar;
        end
    end
    
    out = temp;
end

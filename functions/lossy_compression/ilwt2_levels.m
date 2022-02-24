function out = ilwt2_levels(wimg,wname, levels)

maxlvl = floor(log2(min(size(wimg))));

switch nargin
    case 1
        wname = 'haar';
        levels = maxlvl;
    case 2
        levels = maxlvl;
    case 3
        % do nothing
    otherwise
        warning("Wrong number of arguments");
end

lshaar = liftwave(wname);
maxlvl = floor(log2(min(size(wimg))));

assert( maxlvl >= levels,"Too many levels");


% Add a primal ELS to the lifting scheme.


els = {'p',[-0.125 0.125],0};
lsnew = addlift(lshaar,els);

% inty
lshaarInt = liftwave(wname,'int2int');
lsnewInt = addlift(lshaarInt,els);

orig_size = size(wimg,1);
init_sz = log2(orig_size)-levels;
end_sz = log2(orig_size)-1;

out = wimg;

for ii = init_sz:end_sz
    
    cAintTemp = out(1:2^(ii),1:2^(ii));
    cHintTemp = out(1:2^(ii), 1+2^(ii):2^(ii+1));
    cVintTemp = out(1+2^(ii):2^(ii+1), 1:2^(ii));
    cDintTemp = out(1+2^(ii):2^(ii+1), 1+2^(ii):2^(ii+1));
        
    out(1:2^(ii+1), 1:2^(ii+1))= ...
        ilwt2(cAintTemp,cHintTemp,cVintTemp,cDintTemp, lsnewInt);    
end

function capRGBData =  hs2rgb(hsData, hsBands, refWhite)
% 
%  """Detect specular reflection in a hyperspectral image
% 
% :Example:
% 
%     >>> import idletechs.specular as sp
%     >>> capRGBData = sp.hs2rgb(hsData, hsBands, refWhite)
% 
% :param hsData, input hyperspectral data cube, (nPixels = nSamplesPerLine
% x nLines) x nBands
% :param hsBands, the wavelength to bands map, 1 x nBands
% :param refWhite, reference white for white balancing, [R, G, B]
% :returns: capRGBData, the reconstructed RGB image, (nPixels =
% nSamplesPerLine x nLines) x 3
% :type hsData, numpy.ndarray
% :type hsBands, numpy.ndarray
% :type refWhite, numpy.ndarray
% :rtype capRGBData, numpy.ndarray
% """

assert(ismatrix(hsData));
assert(size(hsData,1)>0);
assert(size(hsData,2)>0);
assert(isnumeric(hsBands));
assert(size(hsBands,1)>0);

nHSBands = size(hsData,2);

if nargin > 2
    assert(isnumeric(refWhite));
    assert(length(refWhite)==3 | length(refWhite)==nHSBands);
end

load('sensitivity.mat','specResp2')

interpSpecResp=zeros(nHSBands,3);

for ch = 1:3
    interpSpecResp(:,ch)=interp1(specResp(1,:),specResp(ch+1,:),hsBands);
end
interpSpecResp(hsBands < specResp(1,1),:) = 0;
interpSpecResp(hsBands > specResp(1,end),:) = 0;


capRGBData = hsData*interpSpecResp;

capRGBData = capRGBData - min(capRGBData);

if nargin > 2
    if length(refWhite)== 3 
            capRGBData = capRGBData./refWhite;
    end

end




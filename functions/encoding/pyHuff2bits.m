%%  Number of bits used in huffman encoding of data 
function [number_of_bits] = pyHuff2bits(data)
% % Input:  data, a vector cointaining real number values to be encoded
% % 
% % Output: number_of_bits, bits used in the encoding with table overhead 

myPyModule = 'huff.py';
pathToHuff = fileparts(which(myPyModule));
if count(py.sys.path,pathToHuff) == 0
    insert(py.sys.path,int32(0),pathToHuff);
end
pyOut = py.huff.huffify(data);

number_of_bits = pyOut.double;
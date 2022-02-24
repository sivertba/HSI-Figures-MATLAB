function [spatialized_scores] = spatialize_scores(scores, cube_size)
% % Input:  spatialized_scores, A numeric matrix representing the scores.
% % Input:  cube_size, A 1x3 vector representing the dimensions of the 
% %         image cube originally
% % 
% % Output: scores, A numeric matrix representing the scores in normally.

assert(ismatrix(scores),"Wrong number of dimensions in scores variable");
assert(length(cube_size) == 3, "Wrong number of dimensions for image cube");

[s1, s2] = size(scores);
assert(min(s1,s2) == s2, "Scores seems to be transposed");

x = cube_size(1);
y = cube_size(2);

score_cube = reshape(scores,[x, y, s2]);

spatialized_scores = zeros(x,y*s2);
for ii = 1:s2
    spatialized_scores(:,((ii-1)*y+1):ii*y) = ...
        score_cube(:,:,ii);
end


end


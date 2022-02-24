function [scores] = despatialize_scores(spatialized_scores, cube_size, score_size)
% % Input:  spatialized_scores, A numeric matrix representing the scores.
% % Input:  cube_size, A 1x3 vector representing the dimensions of the 
% %         image cube originally
% % Input:  score_size, A 1x2 vector representing the dimensions of the 
% %         scores originally
% % 
% % Output: scores, A numeric matrix representing the scores normally.

assert(ismatrix(spatialized_scores),...
    "Wrong number of dimensions in spatialized scores variable");
assert(length(score_size) == 2,...
    "Wrong number of dimensions for score size");
assert(length(cube_size) == 3,...
    "Wrong number of dimensions for cube size");

cx = cube_size(1);
cy = cube_size(2);

s1 = score_size(1);
s2 = score_size(2);
assert(min(s1,s2) == s2, "Scores size seems to be transposed");

temp_cube = zeros(cx, cy, s2);
for ii = 1:s2
    temp_cube(:,:,ii) = ...
    spatialized_scores(:,((ii-1)*cy+1):ii*cy);
end

scores = reshape(temp_cube,[s1, s2]);

end


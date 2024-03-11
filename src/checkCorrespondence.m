function checkCorrespondence(reference, test, flag)

% HOPLAB function to check two different ROIs (one reference and one
% created ad hoc).
% - casts both images as struct type
% - assert equal .mat and .dim
% - returns error in case of irregularities

% If they're not already, cast images as struct
if ~isstruct(reference), reference = struct(reference);
end
if ~isstruct(test), test = struct(test);
end

% Assert equality
% Quick hack to solve for particularity of 'intersection' method
if strcmp(flag, 'intersection')
    assert(isequal(reference.mat, test.mat), 'The affine matrix of the two images is not the same.');
    assert(isequal(reference.dim, test.dim), 'The "size" of the two images is not the same.');
else
    assert(isequal(reference.mat, test.mat), 'The "mat" of the two images is not the same.');
    assert(isequal(reference.dim, size(test.dat)), 'The "size" of the two images is not the same.');
end


end
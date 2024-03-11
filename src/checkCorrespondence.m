function checkCorrespondence(reference, test)

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
assert(isequal(reference.mat, test.mat), 'The "mat" of the two images is not the same.');
assert(isequal(reference.dim, size(test.dat)), 'The "size" of the two images is not the same.');


end
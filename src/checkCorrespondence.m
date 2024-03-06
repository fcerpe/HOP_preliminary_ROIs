function checkCorrespondance(reference, test)

% HOPLAB function to check two different ROIs (one reference and one
% created ad hoc).
% - casts both images as struct type
% - assert equal .mat and .dim
% - returns error in case of irregularities

% cast images as struct
structReference = struct(reference);
structTest = struct(test); 

% Assert equality
assert(isequal(structReference.mat, structTest.mat), 'The "mat" of the two images is not the same.');
assert(isequal(structReference.dim, size(structTest.dat)), 'The "size" of the two images is not the same.');


end
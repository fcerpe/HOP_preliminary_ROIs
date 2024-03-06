function getSelectedVoxels(roiArray, when)

% Print the number of voxels in a ROI

% 'active' ones
numVoxOnes = sum(roiArray(:) == 1);

% 'non-active' ones
numVoxNotOnes = sum(roiArray(:) == 0);

% total number
numVox = length(roiArray(:));
    
% check that total and sums match
assert(isequal(numVoxOnes+numVoxNotOnes, numVox), 'The non-rebased ROI array is not binary.');

% notify the user
fprintf('INFO: selected voxels ', when, ' rebasing: %s\n', num2str(numVoxOnes))


end
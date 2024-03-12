function [thresholdedImage, numCells] = applyThreshold(tConImg, pThresh, erdf)
    % Helper function to apply threshold and count voxels
    tThresh = spm_invTcdf(1-pThresh, erdf); % Calculate T threshold
    thresholdedImage = tConImg > tThresh; % Apply threshold
    numCells = sum(thresholdedImage(:)); % Count voxels above threshold
end
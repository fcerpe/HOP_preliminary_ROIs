function finalRoi = createIntersectedROI(thresholdedImage, roiImg, tStruct, roiName, contrastName)
    % Creates a new ROI object from the intersection of a thresholded T image and an ROI image.
    %
    % This function performs a logical AND operation between a binary thresholded T image
    % and a binary ROI image, marking voxels included in both images. It then creates a new
    % ROI object with a label derived from the ROI name and a simplified contrast name.
    %
    % Parameters:
    %   thresholdedImage: Binary image where voxels above a certain threshold are marked as 1.
    %   roiImg: The ROI image data.
    %   tStruct: The structure returned by spm_vol, containing the T image's affine matrix.
    %   roiName: The name of the ROI.
    %   contrastName: The name of the contrast used for thresholding.
    %
    % Returns:
    %   finalRoi: The newly created ROI object.
    %
    % Example usage:
    %   finalRoi = createIntersectedROI(thresholdedTImg, roiImgData, tImageStruct, 'ROI_Name', 'Contrast_Name');

    fprintf('STEP: Intersecting thresholded T image with ROI...\n');
    
    % Perform intersection by applying a logical AND operation
    intersectedImage = thresholdedImage & (roiImg > 0.5);
    
    % Optionally, count the number of voxels in the intersection for reporting or further analysis
    numIntersectedVoxels = sum(intersectedImage(:));
    fprintf('Number of intersected voxels: %d\n', numIntersectedVoxels);

    fprintf('STEP: Creating new ROI from intersected data...\n');
    
    % Simplify the contrast name for labeling the ROI
    contrastNameSimple = simplifyContrastName(contrastName);
    
    % Create the new ROI object with specified parameters
    finalRoi = maroi_matrix(struct('dat', intersectedImage, 'mat', tStruct.mat, 'label', strcat([roiName, ' ', contrastNameSimple]), 'binarize', 1, 'roithresh', 1e-10));
    
    fprintf('DONE: New ROI created successfully.\n');
end
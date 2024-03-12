function [roiImg, roiStruct] = loadMaskData(roiImgPath)
    % Loads the ROI image and its voxel data from a specified path.
    %
    % Parameters:
    %   roiImgPath: The file path to the ROI image.
    %
    % Returns:
    %   roiImg: The voxel data of the ROI image.
    %   roiStruct: The structure returned by spm_vol, containing image volume information.
    %
    % Example usage:
    %   [roiImg, roiStruct] = loadROIImageData('/path/to/roi/image.nii');
    %
    % This function logs the process of loading ROI data, loads the ROI image,
    % and reads its voxel data, ensuring to handle any errors gracefully.

    fprintf('Loading ROI data from: %s\n', roiImgPath);
    
    try
        % Load the ROI image structure
        roiStruct = spm_vol(roiImgPath);
        
        % Read the voxel data from the ROI image
        roiImg = spm_read_vols(roiStruct);
        
        fprintf('ROI data successfully loaded.\n');
    catch ME
        error('Failed to load ROI data: %s', ME.message);
    end
end

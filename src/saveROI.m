function saveROI(finalRoi, roiImgPath, pThresh, taskName, contrastName, outDir)
    % Saves the ROI to the specified output directory in both NIfTI and .mat formats.
    %
    % This function constructs the output file name using the base ROI name, task name,
    % a simplified contrast name, and the p-value threshold used. It saves the ROI object
    % as an image and as a MATLAB file, logging each step.
    %
    % Parameters:
    %   finalRoi: The ROI object to be saved.
    %   roiImgPath: The file path of the original ROI, used to extract the base file name.
    %   pThresh: The p-value threshold used, for inclusion in the file name.
    %   taskName: The name of the task associated with the ROI.
    %   contrastName: The original contrast name, to be simplified for file naming.
    %   outDir: The directory where the ROI files should be saved.
    %
    % Example usage:
    %   saveROI(finalRoi, '/path/to/original/roi.nii', 0.001, 'task1', 'Contrast 1', '/path/to/output/dir');

    fprintf('STEP: Saving ROI...\n');
    
    % Extract the base file name of the original ROI for use in the output file names
    [~, roiFileName, ~] = fileparts(roiImgPath);

    % Simplify the contrast name for use in the file name
    contrastNameSimple = simplifyContrastName(contrastName);

    % Prepare the p-value string for inclusion in the file name
    pThreshStringSplit = split(string(pThresh), '.');
    pThreshString = pThreshStringSplit(2); % Extracting decimal part for file naming

    % Construct the output file name
    outFileName = strcat(roiFileName, '_task-', taskName, '_contrast-', contrastNameSimple, '_p', pThreshString);

    % Define full paths for the NIfTI and .mat files
    niiFilePath = char(fullfile(outDir, strcat(outFileName, '.nii')));
    matFilePath = char(fullfile(outDir, strcat(outFileName, '.mat')));

    % Save the ROI as an image file
    try
        save_as_image(finalRoi, niiFilePath);
        fprintf('ROI image saved to: %s\n', niiFilePath);
    catch ME
        warning('Failed to save ROI image: %s', ME.message);
    end

    % Save the ROI as a MATLAB .mat file
    try
        saveroi(finalRoi, matFilePath);
        fprintf('ROI MATLAB file saved to: %s\n', matFilePath);
    catch ME
        warning('Failed to save ROI MATLAB file: %s', ME.message);
    end
end
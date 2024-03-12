function [thresholdedImage, finalPThresh] = calculateAndApplyThreshold(tConImg, model, initialPThresh, adjustmentPThreshs, minVoxels)
    % Calculates and applies a threshold to T images based on uncorrected p-values,
    % adjusting the threshold if necessary to meet a minimum voxel count.
    %
    % Parameters:
    %   tConImg: The voxel data of the T contrast image.
    %   model: The loaded SPM model object for obtaining error degrees of freedom.
    %   initialPThresh (optional): The initial p-value threshold. Default is 0.001.
    %   adjustmentPThreshs (optional): Array of p-value thresholds to try if initial thresholding is insufficient. Default is [0.01, 0.05].
    %   minVoxels (optional): Minimum number of voxels required after thresholding. Default is 25.
    %
    % Returns:
    %   thresholdedImage: The T image after applying the final threshold.
    %   finalPThresh: The final p-value threshold used.
    
    % Set default values if not provided
    if nargin < 3 || isempty(initialPThresh)
        initialPThresh = 0.001;
    end
    if nargin < 4 || isempty(adjustmentPThreshs)
        adjustmentPThreshs = [0.01, 0.05];
    end
    if nargin < 5 || isempty(minVoxels)
        minVoxels = 25;
    end
    
    erdf = error_df(model);
    finalPThresh = initialPThresh;
    i = 1; % Start with the first adjustment threshold if needed
    
    % Initially apply threshold
    [thresholdedImage, numCells] = applyThreshold(tConImg, finalPThresh, erdf);
    
    % Adjust threshold if necessary
    while numCells < minVoxels && i <= length(adjustmentPThreshs)
        finalPThresh = adjustmentPThreshs(i);
        [thresholdedImage, numCells] = applyThreshold(tConImg, finalPThresh, erdf);
        i = i + 1; % Move to the next threshold for adjustment
    end
    
    if numCells < minVoxels
        warning('Final voxel count below %d even after adjusting thresholds. Final count: %d', minVoxels, numCells);
    else
        fprintf('Thresholding complete. Final p-value threshold: %f, Voxel count: %d\n', finalPThresh, numCells);
    end
end
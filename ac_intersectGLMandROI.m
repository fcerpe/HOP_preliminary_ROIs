function hop_createROI_intersection(opt, m)

% FUNCTION TO CREATE ROIs FROM AN IMAGE (MASK, SPHERE, PARCEL) AND A CONTRAST ACTIVATION MAP
%
% This script extracts specific regions of interest (ROIs) from an existing 
% ROI image and a statistical contrast activation map previously estimated 
% in SPM. It applies a threshold to the contrast map, masks these voxels to
% the image provided and ensures the final ROI contains a minimum number 
% of voxels. 
% In case the resulting ROI is lower than a specified size, it will relax 
% the threshold until the desired size or a threshold of p < 0.05. 
% Lastly, saves the ROI both in NIfTI and matlab formats.
%
% The script takes the following parameters:
%
%     opt: general options about the folder and files of the project.
%          information about specific paramters can be found in 'hop_opton'
%
%     m: a struct that contains all the specific parameters of this method
%           method: the name of the method requested (e.g. 'atlas')
%           subjects: the list of subjects on which to perform the
%                     extraction. The function expects a valid contrast to
%                     be present (or specified)
%           area: the name of the area that ill be extracted
%           roiPath: the path of the ROI that will be used as mask
%           tmapPath: the path(s) to the contrasted thresholds  
%           task: the task used to find the SPM.mat file containing GLM anlaysis
%           contrast: the name of the contrast to use to extract the
%                     specific t-map
%           targetNbVoxels: the desired minimum number of voxels of the ROI 
%
%
% The script performs the following steps:
%
%     - iterates over each subject
%     - creates an output folder that specifies method and its details, if it does not already exist
%     - loads the sepcified ROI
%     - loads the specified SPM model and contrast
%     - if necessary, resamples the ROIs to match the model space
%     - until the minimum number of voxels is reached, overlap mask and
%       thresholded contrast
%     - binarizes the resulting ROI
%     - saves the merged ROI as a NIfTI file in the output folder
%     - copies the script file to the output folder for replicability.
%
%
% Dependencies: MarsBaR, SPM
%
% Author: Andrea Costantino
% Date: 7 July 2023
%  
% Edited by: Filippo Cerpelloni
% Date: March 2024


% Set parameters
% Define root directories for data, GLM results, MarsBaR toolbox, and ROI paths
derivativesDir = '/data/projects/chess/data/BIDS/derivatives';
GLMroot = fullfile(derivativesDir,'SPM_HPM6_GS1_FD_noHPF/GLM/');
marsabPath = fullfile(derivativesDir, 'marsbar');
roisRoot = fullfile(marsabPath, 'rois-noloc');
outRoot = fullfile(marsabPath, 'rois+loc');

% Setting up ROIs
% FFA - Faces > Objects
roisStruct(1).roiName = 'FFA';
roisStruct(1).roiImgPath = fullfile(opt.dir.output,'method-sphere_radius-10/ROI-FFA_radius-10_space-MNI_resampled-to-sub_binary.nii');
roisStruct(1).taskName = 'loc1';
roisStruct(1).contrastName = 'Faces > Objects';
roisStruct(1).outFolder = 'radius_10+loc';

% LOC - Objects > Scrambled
roisStruct(2).roiName = 'LOC';
roisStruct(2).roiImgPath = fullfile(roisRoot,'radius_10/ROI-LOC_radius-10_space-MNI_resampled-to-sub_binary.nii');
roisStruct(2).taskName = 'loc1';
roisStruct(2).contrastName = 'Objects > Scrambled';
roisStruct(2).outFolder = 'radius_10+loc';

% PPA - Scenes > Objects
roisStruct(3).roiName = 'PPA';
roisStruct(3).roiImgPath = fullfile(roisRoot,'radius_10/ROI-PPA_radius-10_space-MNI_resampled-to-sub_binary.nii');
roisStruct(3).taskName = 'loc1';
roisStruct(3).contrastName = 'Scenes > Objects';
roisStruct(3).outFolder = 'radius_10+loc';

% Find subjects folders
% Define the list of subjects to be processed; use '*' to select all subjects
selectedSubjectsList = '*';

% Get subjects folders from subjects list
subPaths = findSubjectsFolders(GLMroot, selectedSubjectsList);


% Iterate over each subject
for iSub = 1:length(m.subjects)

    % Get subject name
    subName = subPaths(iSub).name;
    fprintf('\n###### Processing %s ######\n', subName);

    % Iterate for each area requested
    for iROI = 1:length(m.roisToCreate)

        % Create output folder
        % skip method details, they will be extracted from the ROI's filename
        outputFolder = createOutputFolder(opt, m, [], subName);

        % Check whether the selected subject has the selected task, otherwise skip
        subPath = fullfile(subPaths(iSub).folder, subPaths(iSub).name); 
        subHasTask = checkSubjectHasTask(subPath, roisStruct(iROI).taskName, subName); % Check if the selected task is present in the sub directory
        if ~subHasTask % If the task name is not found in any of the filename in the subject's directory
            continue; % Skip the rest of the loop and move to the next iteration
        end

        fprintf('\n--- Processing ROI: %s -  %s ---\n',roisStruct(iROI).roiName , roisStruct(iROI).contrastName);

        % Load ROI img data
        [roiImg, roiStruct] = loadROIImageData(roisStruct(iROI).roiImgPath);

        % Load SPM model file
        GLMdir = fullfile(GLMroot, subName, roisStruct(iROI).taskName); % Get SPM folder for this subject and task
        [model, loadSuccess] = loadSPMModel(GLMdir, subName, roisStruct(iROI).taskName); % Import the SPM model file
        if ~loadSuccess % If the SPM model is not found
            continue; % Skip the rest of the loop and move to the next iteration
        end
        
        % Load T constrast and apply threshold (voxel with p < .001)
        [tConImg, tStruct] = loadTContrastImage(model, roisStruct(iROI).contrastName, GLMdir);
        [thresholdedImage, finalPThresh] = calculateAndApplyThreshold(tConImg, model);

        % Check that ROI and t map dimensionality and space are the same
        assert(isequal(roiStruct.mat, tStruct.mat), 'The affine matrix of the two images is not the same.');
        assert(isequal(roiStruct.dim, tStruct.dim), 'The size of the two images is not the same.');

        % Intersect thresholded T map and ROI image
        finalRoi = createIntersectedROI(thresholdedImage, roiImg, tStruct, roisStruct(iROI).roiName, roisStruct(iROI).contrastName);

        % Save the resulting ROI to files(s)
        saveROI(finalRoi, roisStruct(iROI).roiImgPath, finalPThresh, roisStruct(iROI).taskName, roisStruct(iROI).contrastName, outputFolder);
    end

    % Save script
    saveScriptForReplicability(outputFolder)

end

end 




%% Helper functions
function saveScriptForReplicability(outDir)
    % Copies the currently running script to the specified output directory for replicability.
    %
    % This function identifies the full path of the currently executing script and copies it
    % to the given output directory. This process aids in ensuring that analyses can be
    % replicated or reviewed in the future with the exact code version used.
    %
    % Parameters:
    %   outDir: The directory where the script should be copied for future reference.
    %
    % Example usage:
    %   saveScriptForReplicability('/path/to/output/dir');

    fprintf('Saving current script for replicability...\n');
    
    % Get the full path of the currently executing script
    fileNameAndLocation = mfilename('fullpath');
    
    % Extract the directory, file name, and extension of the current script
    [path, filename, ~] = fileparts(fileNameAndLocation);
    
    % Construct the output file name and location
    outputFileNameAndLocation = fullfile(outDir, strcat(filename, '.m'));
    
    % Define the current script's file path
    currentfile = strcat(fileNameAndLocation, '.m');
    
    % Copy the script file to the output directory
    try
        copyfile(currentfile, outputFileNameAndLocation);
        fprintf('Script copied to output folder: %s\n', outputFileNameAndLocation);
    catch ME
        warning('Failed to copy script to output folder: %s', ME.message);
    end
end

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

function contrastNameSimple = simplifyContrastName(contrastName)
    % Simplifies the contrast name by removing non-alphabetic characters and parentheses.
    %
    % Parameters:
    %   contrastName: The original contrast name as a string.
    %
    % Returns:
    %   contrastNameSimple: The simplified contrast name, with only lowercase alphabetic characters.
    %
    % Example usage:
    %   contrastNameSimple = simplifyContrastName('Contrast 1 (Session 1)');
    
    contrastNameSimple = regexprep(lower(contrastName), '[^a-z]+|\([^)]*\)', '');
end

function hasTask = checkSubjectHasTask(subPath, taskName, subName)
    % Checks if the specified subject directory contains files related to the specified task.
    %
    % Parameters:
    %   subPath: The full path to the subject's directory.
    %   taskName: The name of the task to check for.
    %   subName: The name of the subject being checked.
    %
    % Returns:
    %   hasTask: A boolean indicating whether the task is present for the subject.
    %
    % Example usage:
    %   hasTask = checkSubjectHasTask('/path/to/subject/directory', 'task_name', 'subject_name');
    %
    % This function searches the subject's directory for files containing the task name,
    % and provides a warning if the task is not found, suggesting the iteration should be skipped.

    fprintf('Checking for task %s in subject %s directory...\n', taskName, subName);
    
    % Construct the full path and list files
    files = dir(subPath);
    fileNames = {files.name};
    
    % Check for the presence of the task name in any file names
    hasTask = any(contains(fileNames, taskName));
    
    if ~hasTask
        % If the task is not found, issue a warning
        warning('Task %s not found for %s in %s. Skipping...', taskName, subName, subPath);
    else
        fprintf('Task %s found for subject %s.\n', taskName, subName);
    end
end

function [tConImg, tStruct] = loadTContrastImage(model, contrastName, GLMdir)
    % Loads T contrast image voxel data for a given contrast name and GLM directory.
    %
    % Parameters:
    %   model: The loaded SPM model object for the current subject and task.
    %   contrastName: The name of the contrast to load.
    %   GLMdir: Directory containing the GLM results for the subject.
    %
    % Returns:
    %   tConImg: The voxel data of the T contrast image.
    %   tStruct: The structure of the T contrast image.
    %
    % Example usage:
    %   [tConImg, tStruct] = loadTContrastImage(model, 'contrastName', '/path/to/GLMdir');
    %
    % This function first attempts to find the specified contrast within the model.
    % If not found, it appends ' - All Sessions' to the name and retries.
    % It then checks for the existence of the T contrast image file and loads its voxel data.

    try
        fprintf('Attempting to load T contrast image for contrast: %s\n', contrastName);
        
        % Attempt to find the specified contrast by name
        t_con = get_contrast_by_name(model, contrastName);
        if isempty(t_con)
            % If not found, try appending ' - All Sessions' and search again
            fprintf('Contrast %s not found, trying with " - All Sessions" suffix.\n', contrastName);
            contrastName = strcat(contrastName, ' - All Sessions');
            t_con = get_contrast_by_name(model, contrastName);
            if isempty(t_con)
                error('ContrastNotFound', 'Cannot find the contrast %s in the design; has it been estimated?', contrastName);
            end
        end

        % Construct the full path to the T contrast image file
        tConFname = fullfile(GLMdir, t_con.Vspm.fname);

        % Verify the existence of the T contrast image file
        if ~exist(tConFname, 'file')
            error('FileNotFound', 'Cannot find T image %s; has it been estimated?', tConFname);
        else
            fprintf('T contrast image found: %s\n', tConFname);
        end

        % Load the voxel data from the T contrast image
        fprintf('Loading voxel data from T contrast image...\n');
        tStruct = spm_vol(tConFname);
        tConImg = spm_read_vols(tStruct);
        fprintf('Voxel data loaded successfully.\n');

    catch ME
        switch ME.identifier
            case 'ContrastNotFound'
                fprintf('Error: %s\n', ME.message);
            case 'FileNotFound'
                fprintf('Error: %s\n', ME.message);
            otherwise
                fprintf('An unexpected error occurred: %s\n', ME.message);
        end
        tConImg = []; % Return empty in case of error
        tConFname = '';
    end
end

function [model, loadSuccess] = loadSPMModel(GLMdir, subName, taskName)
    % Attempts to load the SPM model for a given subject and task.
    %
    % Parameters:
    %   GLMdir: The directory containing the GLM results for the subject.
    %   subName: The name of the current subject being processed.
    %   taskName: The name of the task for which the model is being loaded.
    %
    % Returns:
    %   model: The loaded SPM model object, or empty if loading failed.
    %   loadSuccess: A boolean indicating whether the model was successfully loaded.
    %
    % Example usage:
    %   [model, loadSuccess] = loadSPMModel('/path/to/GLMdir', 'subject1', 'task1');
    %
    % This function tries to load the SPM.mat file and handles any errors that occur,
    % logging appropriate messages and returning a status indicator.

    fprintf('Attempting to load SPM model for subject %s, task %s...\n', subName, taskName);
    model = []; % Initialize model as empty
    loadSuccess = false; % Initialize success status as false

    try
        % Construct the path to the SPM.mat file
        spmPath = fullfile(GLMdir, 'SPM.mat');

        % Attempt to load the SPM model
        model = mardo(spmPath);
        
        % If successful, set the success status to true
        loadSuccess = true;
        fprintf('SPM model loaded successfully.\n');

    catch ME
        % Handle errors that occur during model loading
        fprintf('WARNING: Error loading SPM model for %s, task %s: %s\n', subName, taskName, ME.message);
        fprintf('Skipping to the next iteration.\n');
        % No need to set model or loadSuccess as they are already initialized to their failure states
    end

    % Return the model (empty if failed) and the success status
    return;
end

function [roiImg, roiStruct] = loadROIImageData(roiImgPath)
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

function [thresholdedImage, numCells] = applyThreshold(tConImg, pThresh, erdf)
    % Helper function to apply threshold and count voxels
    tThresh = spm_invTcdf(1-pThresh, erdf); % Calculate T threshold
    thresholdedImage = tConImg > tThresh; % Apply threshold
    numCells = sum(thresholdedImage(:)); % Count voxels above threshold
end

function [filteredFolderStructure] = findSubjectsFolders(fmriprepRoot, selectedSubjectsList, excludedSubjectsList)
% FINDSUBJECTSFOLDERS Locate subject folders based on a list or wildcard.
%
% USAGE:
% sub_paths = findSubjectsFolders(fmriprepRoot, selectedSubjectsList)
%
% INPUTS:
% fmriprepRoot          - The root directory where 'sub-*' folders are located.
%
% selectedSubjectsList  - Can be one of two things:
%                         1) A list of integers, each representing a subject ID.
%                            For example, [7,9] would search for folders 'sub-07' 
%                            and 'sub-09' respectively.
%                         2) A single character string '*'. In this case, the function
%                            will return all folders starting with 'sub-*'.
%
% OUTPUTS:
% sub_paths             - A structure array corresponding to the found directories.
%                         Each structure has fields: 'name', 'folder', 'date', 
%                         'bytes', 'isdir', and 'datenum'.
%
% EXAMPLES:
% 1) To fetch directories for specific subjects:
%    sub_paths = findSubjectsFolders('/path/to/fmriprepRoot', [7,9]);
%
% 2) To fetch all directories starting with 'sub-*':
%    sub_paths = findSubjectsFolders('/path/to/fmriprepRoot', '*');
%
% NOTE:
% If a subject ID from the list does not match any directory, a warning is issued.

% Start by fetching all directories with the 'sub-*' pattern.
sub_paths = dir(fullfile(fmriprepRoot, 'sub-*'));
sub_paths = sub_paths([sub_paths.isdir]); % Keep only directories.

% Check the type of selectedSubjectsList
if isnumeric(selectedSubjectsList(1))
    % Case 1: selectedSubjectsList is a list of integers.

    % Convert each integer in the list to a string of the form 'sub-XX'.
    subIDs = cellfun(@(x) sprintf('sub-%02d', x), num2cell(selectedSubjectsList), 'UniformOutput', false);

    % Filter the sub_paths to keep only those directories matching the subIDs.
    sub_paths = sub_paths(ismember({sub_paths.name}, subIDs));

    % Check and throw warnings for any missing subID.
    foundSubIDs = {sub_paths.name};
    for i = 1:length(subIDs)
        if ~ismember(subIDs{i}, foundSubIDs)
            warning(['The subID ', subIDs{i}, ' was not found in sub_paths.name.']);
        end
    end

elseif ischar(selectedSubjectsList) && strcmp(selectedSubjectsList, '*')
    % Case 2: selectedSubjectsList is '*'. 
    % No further action required as we've already selected all 'sub-*' folders.

else
    % Invalid input.
    error('Invalid format for selectedSubjects. It should be either "*" or a list of integers.');
end

% Only process exclusion if the excludedSubjectsList is provided.
if nargin == 3
    % Create a list of excluded folder names
    excludedNames = cellfun(@(x) sprintf('sub-%02d', x), num2cell(excludedSubjectsList), 'UniformOutput', false);

    % Logical array of folders to exclude
    excludeMask = arrayfun(@(x) ismember(x.name, excludedNames), sub_paths);

    % Filtered structure
    filteredFolderStructure = sub_paths(~excludeMask);
else
    % If no excludedSubjectsList is provided, just return the sub_paths.
    filteredFolderStructure = sub_paths;
end
end

function hop_roi_intersection(opt, m)

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


% Iterate over each subject
for iSub = 1:length(m.subjects)

    % Get subject name
    subName = ['sub-', m.subjects{iSub}];
    fprintf('\n###### Processing %s ######\n', subName);

    % Iterate for each area requested
    for iROI = 1:length(m.roisToCreate)

        % Take the current ROI that is being extracted
        currROI = m.roisToCreate(iROI);

        % Create output folder
        % method detail is the path to the mask, used to get information
        % from it
        outputFolder = createOutputFolder(opt, m, currROI.maskPath, subName);

        % Notify the user
        fprintf('\n--- Processing ROI: %s -  %s ---\n', currROI.area, currROI.contrast);

        % Load mask image and extract volume data
        [maskImg, maskStruct] = loadMaskData(currROI.maskPath);

        % Get the subject's folder
        % Does not check for presence of the task, node, etc.
        % Refer to 'checkSubjectHasTask' function to improve
        subjectFolder = fullfile(opt.dir.stats, subName, ...
                                ['task-', opt.task,'_space-',opt.space,'_FWHM-', num2str(opt.fwhm),'_node-',opt.node]); 

        % Load SPM model 
        [model, loadSuccess] = loadSPMmodel(subjectFolder, subName, currROI.task); % Import the SPM model file

        % Skip this contrast / ROI if the SPM model fails to load
        if ~loadSuccess, continue 
        end
        
        % Load T constrast and apply threshold (voxel with p < .001)
        [tConImg, tStruct] = loadTContrastImage(model, currROI.contrast, subjectFolder);
        [thresholdedImage, finalPThresh] = calculateAndApplyThreshold(tConImg, model, [], [], m.targetNbVoxels);

        % Check that ROI and t map dimensionality and space are the same
        checkCorrespondence(maskStruct, tStruct, 'intersection');

        % Intersect thresholded T map and ROI image
        finalRoi = createIntersectedROI(thresholdedImage, maskImg, tStruct, currROI.area, currROI.contrast);

        % Save the resulting ROI to files(s)
        saveROI(finalRoi, currROI.maskPath, finalPThresh, currROI.task, currROI.contrast, outputFolder);
    end

    % Save script
    % saveScriptForReplicability(outputFolder)

end

end 



%% Helper functions
% most of them moved to /src as standalone functions
% kept those obsolete or not working

% function saveScriptForReplicability(outDir)
%     % Copies the currently running script to the specified output directory for replicability.
%     %
%     % This function identifies the full path of the currently executing script and copies it
%     % to the given output directory. This process aids in ensuring that analyses can be
%     % replicated or reviewed in the future with the exact code version used.
%     %
%     % Parameters:
%     %   outDir: The directory where the script should be copied for future reference.
%     %
%     % Example usage:
%     %   saveScriptForReplicability('/path/to/output/dir');
% 
%     fprintf('Saving current script for replicability...\n');
%     
%     % Get the full path of the currently executing script
%     fileNameAndLocation = mfilename('fullpath');
%     
%     % Extract the directory, file name, and extension of the current script
%     [path, filename, ~] = fileparts(fileNameAndLocation);
%     
%     % Construct the output file name and location
%     outputFileNameAndLocation = fullfile(outDir, strcat(filename, '.m'));
%     
%     % Define the current script's file path
%     currentfile = strcat(fileNameAndLocation, '.m');
%     
%     % Copy the script file to the output directory
%     try
%         copyfile(currentfile, outputFileNameAndLocation);
%         fprintf('Script copied to output folder: %s\n', outputFileNameAndLocation);
%     catch ME
%         warning('Failed to copy script to output folder: %s', ME.message);
%     end
% end

% function hasTask = checkSubjectHasTask(subPath, taskName, subName)
%     % Checks if the specified subject directory contains files related to the specified task.
%     %
%     % Parameters:
%     %   subPath: The full path to the subject's directory.
%     %   taskName: The name of the task to check for.
%     %   subName: The name of the subject being checked.
%     %
%     % Returns:
%     %   hasTask: A boolean indicating whether the task is present for the subject.
%     %
%     % Example usage:
%     %   hasTask = checkSubjectHasTask('/path/to/subject/directory', 'task_name', 'subject_name');
%     %
%     % This function searches the subject's directory for files containing the task name,
%     % and provides a warning if the task is not found, suggesting the iteration should be skipped.
% 
%     fprintf('Checking for task %s in subject %s directory...\n', taskName, subName);
%     
%     % Construct the full path and list files
%     files = dir(subPath);
%     fileNames = {files.name};
%     
%     % Check for the presence of the task name in any file names
%     hasTask = any(contains(fileNames, taskName));
%     
%     if ~hasTask
%         % If the task is not found, issue a warning
%         warning('Task %s not found for %s in %s. Skipping...', taskName, subName, subPath);
%     else
%         fprintf('Task %s found for subject %s.\n', taskName, subName);
%     end
% end











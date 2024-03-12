function hop_roi_sphere(opt, m)

% FUNCTION TO GENERATE SPHERICAL ROIS IN NIFTI FORMAT
%
% This script generates spherical ROIs (Regions of Interest) based on
% specified MNI coordinates, radii, reference images defined in 'hop_option'. 
% It resamples the ROIs to match a provided reference image and, if asked to, 
% merges them if coordinates are provided for both hemispheres. 
% The resulting ROIs are saved as NIfTI files in an out folder.
%
% The script takes the following parameters:
%
%     opt: general options about the folder and files of the project.
%          information about specific paramters can be found in 'hop_opton'
%
%     m: a struct that contains all the specific parameters of this method
%           method: the name of the method requested (e.g. 'sphere')
%           radii: a vector specifying the radii of the spheres for each ROI
%           area: a list of the names of the areas to extract
%           coordsL: a vector with left-hemisphere MNI coordinates for each area
%           coordsR: a vector with right-hemisphere MNI coordinates for each area
%           referencePath: the path of the reference image used for resampling the ROIs
%           mergeRois: a vector of boolean values that indicates, for each
%                      area, whether to merge the hemispheres in a single ROI
%           
%
% The script performs the following steps:
%
%     - loads the reference image for resampling
%     - iterates over each speficied radius
%     - if not present, creates an output folder for each radius
%     - iterates over each ROI requested
%     - creates sphere ROIs for the left and right hemispheres using the specified MNI coordinates
%     - resamples the ROIs to match the reference image space
%     - verifies the compatibility of the resampled ROI and the reference image
%     - stores the resampled ROIs in a cell array
%     - merges the resampled ROIs if coordinates are provided for both hemispheres
%     - saves the merged ROI as a NIfTI file in the output folder
%     - copies the script file to the output folder for replicability.
%
% The script is executed for each radius specified in the radii vector%
% Dependencies: MarsBaR, SPM
%
% Author: Andrea Costantino
% Date: 6 July 2023
% 
% Edited by: Filippo Cerpelloni
% Date: March 2024


% Iterate for each radius requested
for iRad = 1:length(m.radii)

    % Pick the current radius
    radius = m.radii(iRad);

    % Notify the user
    fprintf(['\nExtracting ROIs for ', num2str(radius),'mm radius\n']);

    % Create output folder and pass the radius as argument
    outputFolder = createOutputFolder(opt, m, radius);

    % Create a mars_space object for the reference image, and save a copy
    % in the output folder
    [referenceSpace, ~] = loadReference(opt, m, outputFolder);


    % Loop through each area specified for this method
    % TL;DR
    % - extract ROI
    % - resample it on reference image from the corresponding space
    % - merge hemispheres if required 
    % - save all ROIs in the output folder
    for iRTC = 1:length(m.roisToCreate)

        % Select current ROI to work on
        currROI = m.roisToCreate(iRTC);

        % Notifiy user
        fprintf(['\nSTEP: Processing ROI: ', currROI.area, '\n']);

        % If both hemisphere coordinates have AT LEAST one NaN value, skip this ROI 
        if isempty(currROI.coordsL) && isempty(currROI.coordsR)
            
            % Throw a warning, there may be a problem in the options
            warning('Skipping this area: missing some coordinates. Check options to make sure everyhting is in order');
            continue
        end

        % Initialize empty cell array to store resampled ROIs
        processedRoiList = {};

        % For each set of MNI coordinates (L and R), create sphere ROI
        % and resample to match the reference image
        for hemi = ['L', 'R']

            % Fetch the current hemisphere
            coords = currROI.(strcat('coords', hemi));

            % If it's defined
            if ~any(isnan(coords)) && ~isempty(coords)

                % Notify the user
                fprintf(['  STEP: Creating sphere ROI for ', hemi, ' hemisphere\n']);

                % Create the sphere ROI
                sphere_roi = maroi_sphere(struct('centre', coords, 'radius', radius, ...
                                          'label', strcat('hemi-', hemi, '_label-', currROI.area), ...
                                          'binarize', 1, 'roithresh', 0.5));

                % Resample the ROI
                resampledRoi = maroi_matrix(sphere_roi, referenceSpace);

                % Check mat and sizes of the resampled ROI against the reference image
                checkCorrespondance(referenceSpace, resampledRoi, 'sphere');

                % Store the resampled ROI in the cell array
                processedRoiList{end + 1} = resampledRoi;

                fprintf(['DONE: Resampled ROI created for ', hemi, ' hemisphere\n']);
            end
        end

        % If required, merge the resampled ROIs
        if currROI.mergeRois && length(processedRoiList) > 1
            fprintf('STEP: Merging resampled ROIs\n');
        
            % Make empty array with same dim as refImg
            mergedRoiMat = zeros(referenceSpace.dim);

            % Loop through and merge remaining ROIs
            for iROI = 1:length(processedRoiList)

                % Merge current ROI with the final ROI
                currentRoiStruct = struct(processedRoiList{iROI}).dat;
                mergedRoiMat = currentRoiStruct > 0.5  | mergedRoiMat > 0.5;
            end
                
            % Make a new ROI with same affine and dimension of the refImg
            mergedRoi = maroi_matrix(struct('dat', mergedRoiMat, 'mat', referenceImage.mat, 'label', ...
                                    strcat('hemi-B_label-', m.area{iRTC}), ...
                                    'binarize', 1, 'roithresh', 0.5 ));

            fprintf('DONE: ROIs merged successfully!\n');

            % Quick hack: rename 'hemi' to [B]ilateral, to be used in filename
            hemi = 'B';

        else
            % If there's only one resampled ROI, simply use it as the final ROI
            fprintf('SKIP: Only one resampled image provided. Skipping merging hemispheres...\n');
            mergedRoi = processedRoiList{1};
        end

        % Verify that the final (merged, resampled) ROI and the reference image are in the same space
        checkCorrespondence(referenceSpace, mergedRoi, 'sphere');
        

        % Save the merged ROI
        fprintf('STEP: Saving ROIs\n');
        filename = fullfile(outputFolder, ['hemi-', hemi, '_space-', opt.space, ...
                                           '_method-', m.method, '_radius-', num2str(radius), ...
                                           '_label-', currROI.area, '_desc-resampled_mask']);

        save_as_image(mergedRoi, [filename, '.nii']);
        saveroi(mergedRoi, [filename, '.mat']);

        fprintf('DONE: ROI saved as: %s\n', filename);
    end

    fprintf(['\nExtracted all the ROIs for ', num2str(radius),'mm radius\n']);

    % Save script - ON HOLD - DOES NOT WORK
    % Copy this script to the output folder for replicability
%     scriptPath = mfilename('fullpath');
%     [path, filename, ext] = fileparts(scriptPath);
%     scriptSaveFilenmame = fullfile(opt.dir.jobs, strcat(filename, '.m'));
%     currentfile = strcat(fileNameAndLocation, '.m');
%     copyfile(currentfile, outputFileNameAndLocation);
%     fprintf('Script copied to output folder\n');

end

end
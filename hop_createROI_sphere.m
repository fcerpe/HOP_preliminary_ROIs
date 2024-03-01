function hop_createROI_sphere(opt, m)

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
%     - iterates over 
%     - if not present, creates an output folder for each radius
%     - iterates over each ROI in m
%     -  and creates a mars_space object for resampling
%     - creates sphere ROIs for the left and right hemispheres using the specified MNI coordinates
%     - resamples the ROIs to match the reference image space
%     - verifies the compatibility of the resampled ROI and the reference image
%     - stores the resampled ROIs in a cell array
%     - merges the resampled ROIs if coordinates are provided for both hemispheres
%     - saves the merged ROI as a NIfTI file in the output folder
%     - copies the script file to the output folder for replicability.
%
% The script is executed for each radius specified in the radii vector.%
% Dependencies: MarsBaR, SPM
%
% Author: Andrea Costantino
% Date: 6 July 2023
% 
% Edited by: Filippo Cerpelloni
% Date: March 2024


% Load reference image provided by method, to be used for resampling
referenceImg = spm_vol(m.referencePath);
referenceImg = referenceImg(1, 1);

% Create a mars_space object for the reference image
referenceSpace = mars_space(referenceImg);


for iRad = 1:length(m.radii)

    % Pick the current radius
    radius = m.radii(iRad);

    % Notify the user
    fprintf(['\nExtracting ROIs for ', num2str(radius),'mm radius\n']);

    % If required, create output folder. In intersectGLMandROI there is a
    % subfuction for this
    % add information about method used and method-specific parameters
    outPath = fullfile(opt.dir.output, ['method-', m.method '_radius-', num2str(radius)]);
    mkdir(outPath)

    fprintf(['ROIs will be saved in: ''', outPath, '''\n']);


    % Loop through each area specified for this method
    % TL;DR
    % - extract ROI
    % - resample it on reference image from the corresponding space
    % - merge hemispheres if required 
    % - save all ROIs in the output folder
    for iArea = 1:length(m.area)

        % Notifiy user
        fprintf(['\nSTEP: Processing ROI: ', m.area{iArea}, '\n']);

        % If both hemisphere coordinates have AT LEAST one NaN value, skip this ROI 
        if any(isnan(m.coordsL(iArea))) && any(isnan(m.coordsR(iArea)))
            
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
            coords = m.(strcat('coords', hemi));

            % If it's defined
            if ~any(isnan(coords(iArea)))

                % Notify the user
                fprintf(['  STEP: Creating sphere ROI for ', hemi, ' hemisphere\n']);

                % Create the sphere ROI
                sphere_roi = maroi_sphere(struct('centre', coords, 'radius', radius, ...
                                          'label', strcat('hemi-', hemi, '_label-', m.area{iArea}), ...
                                          'binarize', 1, 'roithresh', 0.5));

                % Resample the ROI
                resampledRoi = maroi_matrix(sphere_roi, referenceSpace);

                % Check mat and sizes of the resampled ROI against the reference image
                processedRoi = struct(resampledRoi);
                assert(isequal(referenceImg.mat, processedRoi.mat), 'The "mat" of the two images is not the same.');
                assert(isequal(referenceImg.dim, size(processedRoi.dat)), 'Sizes are not the same');

                % Store the resampled ROI in the cell array
                processedRoiList{end + 1} = resampledRoi;

                fprintf(['DONE: Resampled ROI created for ', hemi, ' hemisphere\n']);
            end
        end

        % If required, merge the resampled ROIs
        if m.mergeRois(iArea) && length(processedRoiList) > 1
            fprintf('STEP: Merging resampled ROIs\n');
        
            % Make empty array with same dim as refImg
            mergedRoiMat = zeros(referenceImg.dim);

            % Loop through and merge remaining ROIs
            for iROI = 1:length(processedRoiList)

                % Merge current ROI with the final ROI
                currentRoiStruct = struct(processedRoiList{iROI}).dat;
                mergedRoiMat = currentRoiStruct > 0.5  | mergedRoiMat > 0.5;
            end
                
            % Make a new ROI with same affine and dimension of the refImg
            mergedRoi = maroi_matrix(struct('dat', mergedRoiMat, 'mat', referenceImg.mat, 'label', ...
                                    strcat('hemi-B_label-', m.area{iArea}), ...
                                    'binarize', 1, 'roithresh', 0.5 ));

            fprintf('DONE: ROIs merged successfully!\n');

            % Quick hack: modify 'hemi' to [B]ilateral, to be used in
            % filename
            hemi = 'B';

        else
            % If there's only one resampled ROI, simply use it as the final ROI
            fprintf('SKIP: Only one resampled image provided. Skipping merging hemispheres...\n');
            mergedRoi = processedRoiList{1};
        end

        % Verify that the final (merged, resampled) ROI and the reference image are in the same space
        mergedRoiStruct = struct(mergedRoi);
        assert(isequal(referenceImg.mat, mergedRoiStruct.mat), 'The "mat" of the two images is not the same.');
        assert(isequal(referenceImg.dim, size(mergedRoiStruct.dat)), 'The "size" of the two images is not the same.');


        % Save the merged ROI
        fprintf('STEP: Saving ROIs\n');
        filename = fullfile(outPath, ['hemi-', hemi, '_space-', opt.space, ...
                                      '_method-sphere_radius-10_label-', m.area{iArea}, '_desc-resampled-to-sub_binary_mask']);

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
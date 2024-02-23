function ac_makeSphereROI(opt, m)

% SCRIPT FOR GENERATING SPHERICAL ROI IN NIFTI FORMAT
%
% This script generates spherical ROIs (Regions of Interest) based on
% specified MNI coordinates and radii. It resamples the ROIs to match a
% reference image and merges them if coordinates are provided for both
% hemispheres. The resulting ROIs are saved as NIfTI files in an out folder.
%
% The script takes the following parameters:
%
%     radii: A vector specifying the radii of the spheres for each ROI.
%     roisStruct: A structure array containing information about each ROI, including the ROI name and MNI coordinates for the left and right hemispheres.
%               e.g., % FFA ROI
%                     roisStruct(1).roiName = 'FFA';
%                     roisStruct(1).mni_coordinates_l = [-38, -58, -14];
%                     roisStruct(1).mni_coordinates_r = [40, -55, -12];
%
%     refImagePath: The path of the reference image used for resampling the ROIs.
%     outRoot: The root folder where the output will be saved.
%
% The script performs the following steps:
%
%     - Creates an output folder for each radius.
%     - Iterates over each ROI in the roisStruct.
%     - Loads the reference image and creates a mars_space object for resampling.
%     - Creates sphere ROIs for the left and right hemispheres using the specified MNI coordinates.
%     - Resamples the ROIs to match the reference image space.
%     - Verifies the compatibility of the resampled ROI and the reference image.
%     - Stores the resampled ROIs in a cell array.
%     - Merges the resampled ROIs if coordinates are provided for both hemispheres.
%     - Saves the merged ROI as a NIfTI file in the output folder.
%     - Copies the script file to the output folder for replicability.
%
% The script is executed for each radius specified in the radii vector.%
% Dependencies: MarsBaR, SPM
%
% Author: Andrea Costantino
% Date: 6 July 2023


% Parameters
% Define the ROIs
radii = m.radii;

roisStruct = struct(); % Initialize the structure

% FFA ROI
roiName = 'FFA';
mni_coordinates_l = [-38, -58, -14];
mni_coordinates_r = [40, -55, -12];

% Reference image path (subject's functional image)
refImagePath = '/data/projects/chess/data/BIDS/derivatives/SPM/gunzipped/sub-05/sub-05_task-exp_run-1_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii';      % Path of the reference image


for iR = 1:length(m.radii)
    radius = m.radii(iR);

    fprintf('\n###### Processing Radius %s ######\n', num2str(radius));

    % If required, create output folder
    if opt.saveInSubfolder

        outPath = fullfile(opt.dir.output, ['method-', m.method '_radius-', num2str(radius)]);
        mkdir(outPath)
        fprintf('DONE: Created output folder: %s\n', outPath);
    end

    % Make ROI, resample, and merge if two hemispheres are provided
    % Loop through each ROI in the structure
    for iROI = 1:length(roisStruct)


        fprintf('\n--- Processing ROI: %s -  %s ---\n',roiName , num2str(radius));

        % Skip this iteration if both coordinates are empty
        if isempty(mni_coordinates_l) && isempty(mni_coordinates_r)
            fprintf('SKIP: Both hemisphere coordinates are empty. Skipping...\n');
            continue;
        end

        % Load reference image for resampling
        refImgs = spm_vol(refImagePath);
        refImg = refImgs(1, 1);

        % Create a mars_space object for the reference image
        refSpace = mars_space(refImg);

        % Initialize empty cell array to store resampled ROIs
        resampledRoisStruct = {};

        % For each set of MNI coordinates (L and R), create sphere ROI
        % and resample to match the reference image
        for side = ["l", "r"]
            mni_coordinates = roiStruct.(strcat('mni_coordinates_', side));
            if ~isempty(mni_coordinates)
                fprintf('STEP: Creating sphere ROI for %s hemisphere\n', side);

                % Create the sphere ROI
                sphere_roi = maroi_sphere(struct('centre', mni_coordinates, 'radius', radius, 'label', strcat(roiName, '_', side), 'binarize', 1, 'roithresh', 0.5));

                % Resample the ROI
                resampledRoi = maroi_matrix(sphere_roi, refSpace);

                % Verify that the resampled ROI and the reference image are in the same space
                resampledRoiStruct = struct(resampledRoi);
                assert(isequal(refImg.mat, resampledRoiStruct.mat), 'The "mat" of the two images is not the same.');
                assert(isequal(refImg.dim, size(resampledRoiStruct.dat)), 'The "size" of the two images is not the same.');

                % Store the resampled ROI in the cell array
                resampledRoisStruct{end + 1} = resampledRoi;

                fprintf('DONE: Resampled ROI created for %s hemisphere\n', side);
            end
        end

        % Merge the resampled ROIs
        if length(resampledRoisStruct) > 1
            fprintf('STEP: Merging resampled ROIs\n');
        
            % Make empty array with same dim as refImg
            finalRoiMat = zeros(refImg.dim);

            % Loop through and merge remaining ROIs
            for iROI = 1:length(resampledRoisStruct)
                % Merge current ROI with the final ROI
                currentRoiStruct = struct(resampledRoisStruct{iROI}).dat;
                finalRoiMat = currentRoiStruct > 0.5  | finalRoiMat > 0.5;
            end
                
            % Make a new ROI with same affine and dimension of the refImg
            finalRoi = maroi_matrix(struct('dat', finalRoiMat, 'mat', refImg.mat, 'label', roiName, 'binarize', 1, 'roithresh', 0.5 ));

            fprintf('DONE: ROIs merged successfully!\n');
        else
            % If there's only one resampled ROI, simply use it as the final ROI
            fprintf('SKIP: Only one resampled image provided. Skipping merging hemispheres...\n');
            finalRoi = resampledRoisStruct{1};
        end

        % Verify that the final (merged, resampled) ROI and the reference image are in the same space
        finalRoiStruct = struct(finalRoi);
        assert(isequal(refImg.mat, finalRoiStruct.mat), 'The "mat" of the two images is not the same.');
        assert(isequal(refImg.dim, size(finalRoiStruct.dat)), 'The "size" of the two images is not the same.');

        % Save the merged ROI
        fprintf('STEP: Saving merged ROIs\n');
        filePath = fullfile(outPath, strcat('ROI-', roiName, '_','radius-', num2str(radius), '_space-MNI_resampled-to-sub_binary'));

        save_as_image(finalRoi, strcat(filePath, '.nii'));
        saveroi(finalRoi, strcat(filePath, '.mat'));

        fprintf('DONE: Merged ROI saved as: %s\n', filePath);
    end

    fprintf('###### Radius %s processing completed ######\n\n', num2str(radius));

    % Save script
    % Copy this script to the output folder for replicability
    fileNameAndLocation = [mfilename('fullpath')];
    [path, filename, ext] = fileparts(fileNameAndLocation);
    outputFileNameAndLocation = fullfile(outPath, strcat(filename, '.m'));
    currentfile = strcat(fileNameAndLocation, '.m');
    copyfile(currentfile, outputFileNameAndLocation);
    fprintf('Script copied to output folder\n');

end

end
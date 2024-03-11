% SCRIPT FOR EXTRACTING AND RESAMPLING ROIS FROM A BRODMANN ATLAS
%
% This script extracts specific regions of interest (ROIs) from a Brodmann area (BA) atlas and resamples them to match a given reference image. The resulting resampled ROIs are then saved as NIfTI and MATLAB files in the specified output folder.
%
% The script takes the following parameters:
%
% roisStruct: A structure array specifying the ROIs to extract. Each element of the array should have the following fields:
% roiName: A string that specifies the name of the ROI.
% BAs: A vector of integers specifying the Brodmann areas to include in the ROI.
% Example:
% roisStruct(1).roiName = 'DLPFC';
% roisStruct(1).BAs = [9, 46];
%
% baAtlasPath: A string that specifies the path to the Brodmann area atlas in NIfTI format.
%
% refImagePath: A string that specifies the path to the reference image that the extracted ROIs will be resampled to match.
%
% outRoot: A string that specifies the root directory where the output files will be saved.
%
% The script performs the following steps:
%
% - It creates an output directory if it does not already exist.
% - It imports the reference image, decompresses it if necessary, and extracts the data and affine transformation from it.
% - It imports the Brodmann area atlas, decompresses it if necessary, and extracts the voxel data from it.
% - For each ROI specified in roisStruct, it:
% - Extracts the voxels corresponding to the specified Brodmann areas from the atlas.
% - Binarizes the resulting mask.
% - Resamples the binarized mask to match the space of the reference image.
% - Saves the resampled ROI as a NIfTI file and a MATLAB file in the output directory.
% - It saves a copy of the script in the output directory for reproducibility purposes.
%
% Dependencies: MarsBaR, SPM
%
% Author: Andrea Costantino
% Date: 7 July 2023

clear
clc

%% Parameters
% Define the ROIs
roisStruct = struct(); % Initialize the structure

% DLPFC ROI
roisStruct(1).roiName = 'DLPFC';
roisStruct(1).BAs = [9, 46];

% LVC ROI
roisStruct(2).roiName = 'LVC';
roisStruct(2).BAs = [17, 18];

% MRIcro BA atlas path
baAtlasPath = '/usr/share/mricron/templates/brodmann.nii.gz';

% Reference image path (subject's functional image)
refImagePath = '/data/projects/chess/data/BIDS/derivatives/SPM/gunzipped/sub-05/sub-05_task-exp_run-1_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii';      % Path of the reference image

% Set output folder
outRoot = '/data/projects/chess/data/BIDS/derivatives/marsbar';

%% Create output folder
% Set and create the output path
outPath = fullfile(outRoot, strcat('BA'));
if ~exist(outPath, 'dir')
    fprintf('STEP: Creating output folder\n');
    mkdir(outPath)
    fprintf('DONE: Created output folder: %s\n', outPath);
else
    fprintf('SKIP: Output folder already exists: %s\n', outPath);
end

%% Import reference image, gunzip if needed
[refPath, refFilename, refExt] = fileparts(refImagePath);
if strcmp(refExt, '.gz')
    % Gunzip and extract gunzipped file path
    refPathCell = gunzipNiftiFile(refImagePath, outRoot);
    refImagePath = refPathCell{1};

    % Get the filename without the .nii extension
    atlasFilenameSplit = split(refFilename, '.');
    refFilename = atlasFilenameSplit{1};
elseif strcmp(refExt, '.nii')
    % Copy the file from the source folder to the destination folder
    newrefImagePath = fullfile(outPath, strcat(refFilename, '.nii'));
    copyfile(refImagePath, newrefImagePath);
    refImagePath = newrefImagePath;
else
    error('No file ATLAS file found!');
end

% Load reference image for resampling
fprintf('STEP: Getting REF data \n');
refImgsStructLong = spm_vol(refImagePath);
refImgStruct = refImgsStructLong(1,1);
refSpace = mars_space(refImgStruct);
fprintf('DONE: REF data loaded \n');

%% Import BA atlas and gunzip if needed
fprintf('STEP: Getting ATLAS data \n');
[atlasPath, atlasFilename, atlasExt] = fileparts(baAtlasPath);
if strcmp(atlasExt, '.gz')
    % Gunzip and extract gunzipped file path
    baAtlasPathCell = gunzipNiftiFile(baAtlasPath, outRoot);
    baAtlasPath = baAtlasPathCell{1};

    % Get the filename without the .nii extension
    atlasFilenameSplit = split(atlasFilename, '.');
    atlasFilename = atlasFilenameSplit{1};
elseif strcmp(atlasExt, '.nii')
    % Copy the file from the source folder to the destination folder
    newBaAtlasPath = fullfile(outPath, strcat(atlasFilename, '.nii'));
    copyfile(baAtlasPath, newBaAtlasPath);
    baAtlasPath = newBaAtlasPath;
else
    error('No file ATLAS file found!');
end

% Load atlas img and voxel data
fprintf('STEP: Getting ATLAS data \n');
baAtlasStruct = spm_vol(baAtlasPath);
baAtlas = spm_read_vols(baAtlasStruct);

fprintf('DONE: ATLAS data loaded.\n');

%% Make ROI, resample, and merge if two labels are provided
% Loop through each ROI in the structure
for roiID = 1:length(roisStruct)

    % Retrieve the current ROI's parameters
    roiStruct = roisStruct(roiID);
    roiName = roiStruct.roiName;
    roiBAs = roiStruct.BAs;
    BAstring = join(split(num2str(roiBAs),'  '), '+');

    fprintf('\n--- Processing ROI: %s ---\n',roiName);

    % Skip this iteration if both coordinates are empty
    if isempty(roiBAs)
        fprintf('SKIP: Brodmann areas labels are empty. Skipping...\n');
        continue;
    end

    % Initialize empty array to store selected voxels
    roiArray = zeros(baAtlasStruct.dim);

    % For each BA label get the voxels in that area
    for baID = 1:length(roiBAs)
        BA = roiBAs(baID);
        if ~isempty(BA)

            fprintf('STEP: Extracting area %s for ROI %s\n', num2str(BA), roiName);
            % Assign 1 to the final array if the array is 1 or the atlas
            % voxels belong to the selected Brodmann area
            roiArray = roiArray == 1 | baAtlas == BA;

        end
    end

    % Binarize the mask after rebasing
    finalRoiArray = double(roiArray);

    % Get and print number of selected voxels ( == 1)
    numVoxOnes = sum(finalRoiArray(:) == 1);
    numVoxNotOnes = sum(finalRoiArray(:) == 0);
    numVox = length(finalRoiArray(:));
    
    assert(isequal(numVoxOnes+numVoxNotOnes, numVox), 'The non-rebased ROI array is not binary.');

    fprintf('INFO: selected voxels BEFORE rebasing: %s\n', num2str(numVoxOnes))

    % Now that we have the final array, create a maroi object to resample 
    % the final array to our desired space (the subject's space -> MNI)
    roiMaroi = maroi_matrix(struct('dat', finalRoiArray, 'mat', baAtlasStruct.mat, 'label', sprintf('Binarized %s [%s]', roiName, num2str(roiBAs)), 'binarize', 1, 'roithresh', 0.5 ));
    
    % Resample and reslice the ROI
    resampledRoi = maroi_matrix(roiMaroi, refSpace);

    % Verify that the resampled ROI and the reference image are in the same space
    resampledRoiStruct = struct(resampledRoi);

    assert(isequal(refImgStruct.mat, resampledRoiStruct.mat), 'The "mat" of the two images is not the same.');
    assert(isequal(refImgStruct.dim, size(resampledRoiStruct.dat)), 'The "size" of the two images is not the same.');

    % Get and print number of selected voxels after rebasing ( == 1)
    numResampledVoxOnes = sum(resampledRoiStruct.dat(:) == 1);
    numResampledVoxNotOnes = sum(resampledRoiStruct.dat(:) == 0);
    numResampledVox = length(resampledRoiStruct.dat(:));

    assert(isequal(numVoxOnes+numVoxNotOnes, numVox), 'The rebased ROI array is not binary.');

    fprintf('INFO: selected voxels AFTER rebasing: %s\n', num2str(numResampledVoxOnes))

    % Save resampled ROI
    outRoiName = sprintf('ROI-%s_BA-%s_binarized', roiName, BAstring{1});
    outRoiPath = fullfile(outPath, outRoiName);

    save_as_image(resampledRoi, strcat(outRoiPath, '.nii'));
    saveroi(resampledRoi, strcat(outRoiPath, '.mat'));

    fprintf('DONE: Resampled %s ROI saved as: %s\n', roiName, outRoiPath);

    fprintf('Processing ROI %s: COMPLETED!\n\n',roiName);
            
end

%% Save script
% Copy this script to the output folder for replicability
fileNameAndLocation = [mfilename('fullpath')];
[path, filename, ext] = fileparts(fileNameAndLocation);
outputFileNameAndLocation = fullfile(outPath, strcat(filename, '.m'));
currentfile = strcat(fileNameAndLocation, '.m');
copyfile(currentfile, outputFileNameAndLocation);
fprintf('Script copied to output folder\n');


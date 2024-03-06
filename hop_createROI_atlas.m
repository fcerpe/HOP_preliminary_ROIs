function hop_createROI_atlas(opt, m)

% FUNCTION TO EXTRACT ROIS FROM A GIVEN ATLAS
%
% This script extracts specific regions of interest (ROIs) from provided 
% parcels of an atlas and resamples them to match a given reference image. 
% The resulting resampled ROIs are then saved as NIfTI and MATLAB files in 
% the specified output folder.
%
% The script takes the following parameters:
%
%     opt: general options about the folder and files of the project.
%          information about specific paramters can be found in 'hop_opton'
%
%     m: a struct that contains all the specific parameters of this method
%           method: the name of the method requested (e.g. 'atlas')
%           referencePath: the path of the reference image used for resampling the ROIs
%           atlas: the name of the atlas used 
%           mask: an array containing the label(s) of the parcel(s) to
%                 extract
%           roiPath: an array containing the paths of the masks
%
%
% The script performs the following steps:
%
%     - creates an output folder that specifies method and atlas used, if it does not already exist
%     - loads the reference image for resampling
%     - loads the atlas requested, and extracts voxel data 
%     - iterates over each ROI requested
%     - extracts the voxels corresponding to the specified parcel(s)
%     - binarizes the resulting ROI
%     - resamples the ROIs to match the reference image space
%     - verifies the compatibility of the resampled ROI and the reference image
%     - stores the resampled ROIs in a cell array
%     - merges the resampled ROIs if multiple parcels are provided
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

% Prepare inputs and outputs

% Create output folder
% will create a subfolder in output with the format 'mathod-atlas_atlas-[atlas specified]
% If the folder is already present, will not overwrite. Be careful not to
% overwrite files
outputFolder = createOutputFolder(opt, m);

% Import the reference image
% - gunzip if necessary
% - save a copy in the method's output folder
[referenceSpace, referenceStruct] = loadReference(opt, m, outputFolder);

% Import atlas
% - from atlas requested, pick the correct path
% - gunzip if necessary
% - save a copy in the method's output folder
[atlas, atlasStruct] = loadAtlas(m, outputFolder);


% Loop through each ROI requested
% - make ROI 
% - resample
% - merge if two labels are provided
for iRTC = 1:length(m.roisToCreate)

    % Select current ROI to work on
    currROI = m.roisToCreate(iRTC);

    % Notifiy user
    fprintf(['\nSTEP: Processing ROI: ', currROI.area, '\n']);

     
    % If method is custom, there won't be a parcels field in the struct.
    % We can assign '1' assuming it will be a regular mask
    if ~isfield(currROI, 'parcels'), currROI.parcels = 1;
    end
 
    % Extract parcels
    parcelsString = join(split(num2str(currROI.parcels),'  '), '-');

    % Notify the user
    fprintf('\n--- Processing ROI: %s ---\n', currROI.area);

    % Skip this iteration if both coordinates are empty
    if isempty(currROI.parcels)
        fprintf('No parcels provided. Skipping...\n');
        continue
    end

    % Initialize empty array to store selected voxels
    roiArray = zeros(atlasStruct.dim);

    % For each label get the voxels in that area
    for iPar = 1:length(currROI.parcels)
        parcel = currROI.parcels(iPar);
        if ~isempty(parcel)

            fprintf('STEP: Extracting area %s for ROI %s\n', num2str(parcel), currROI.area);

            % Assign 1 to the final array if the array is 1 or the atlas
            % voxels belong to the selected area
            roiArray = roiArray == 1 | atlas == parcel;
        end
    end

    % Binarize the mask after rebasing
    jointRoiArray = double(roiArray);

    % Get and print number of selected voxels 
    getSelectedVoxels(jointRoiArray, 'before');

    % Now that we have the final array, create a maroi object to resample 
    % the final array to our desired space (the subject's space -> MNI)
    marsRoi = maroi_matrix(struct('dat', jointRoiArray, 'mat', atlasStruct.mat, 'label', ...
                                   sprintf('Binarized %s [%s]', currROI.area, num2str(currROI.parcels)), 'binarize', 1, 'roithresh', 0.5 ));
    
    % Resample and reslice the ROI
    resampledRoi = maroi_matrix(marsRoi, referenceSpace);

    % Verify that the resampled ROI and the reference image are in the same space
    checkCorrespondence(referenceSpace, resampledRoi);

    % Get and print number of selected voxels after rebasing ( == 1)
    resampledStruct = struct(resampledRoi);
    getSelectedVoxels(resampledStruct.dat, 'after')


    % Save resampled ROI
    fprintf('STEP: Saving ROIs\n');
    filename = fullfile(outputFolder, ['space-', opt.space, ...
                                       '_method-', m.method, '_atlas-', m.atlas, ...
                                       '_label-', currROI.area, '_parcels-', parcelsString{1}, '_mask']);

    save_as_image(resampledRoi, strcat(filename, '.nii'));
    saveroi(resampledRoi, strcat(filename, '.mat'));

    fprintf('DONE: Resampled %s ROI saved as: %s\n', currROI.area, filename);

    fprintf('Processing ROI %s: COMPLETED!\n\n',currROI.area);
            
end

% Save script
% Copy this script to the output folder for replicability
% fileNameAndLocation = [mfilename('fullpath')];
% [path, filename, ext] = fileparts(fileNameAndLocation);
% outputFileNameAndLocation = fullfile(outPath, strcat(filename, '.m'));
% currentfile = strcat(fileNameAndLocation, '.m');
% copyfile(currentfile, outputFileNameAndLocation);
% fprintf('Script copied to output folder\n');

end
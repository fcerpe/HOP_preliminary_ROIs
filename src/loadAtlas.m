function [atlas, atlasStruct] = loadAtlas(m, outputFolder)

% HOPLAB function to load a given atlas for ROI creation
% - assumes a valid atlas is provided as input
% - if 'custom' is the name of the atlas, loads the .nii file without any
%   checks


% pick the path based on the user's request
switch m.atlas

    case 'Brodmann'
        atlasPath = 'masks/Brodmann/brodmann.nii.gz';

    % Treat the custom atlas as a regular one, but specify the path to
    % be the one provided as parameter
    case 'custom'
        atlasPath = m.roiPath;
end


% Notify the user
fprintf('STEP: Getting ATLAS data \n');

% Extract parts of the filename
[atPath, atFilename, atExt] = fileparts(atlasPath);

newAtlasName = ['reference_atlas-', m.atlas, '.nii'];

% If needed, gunzip NIfTI file 
if strcmp(atlasExt, '.gz')

    % Gunzip and overwrite path with new one
    atlasPath = gunzip(atlasPath, optputFolder);
    atlasPath = atlasPath{1};

    % Get the filename without the .nii extension
    atlasFilenameSplit = split(atFilename, '.');
    atlasFilename = atlasFilenameSplit{1};

elseif strcmp(atlasExt, '.nii')

    % Copy the file from the source folder to the destination folder
    newAtlasPath = fullfile(outputFolder, newAtlasName);
    copyfile(atlasPath, newAtlasPath);
    atlasPath = newAtlasPath;

else
    error('No file ATLAS file found!');
end

% Load atlas img and voxel data

% Notifiy the user
fprintf('STEP: Getting ATLAS data \n');

% Load the NIfTI file with SPM functions, stor both the struct and the space
atlasStruct = spm_vol(atlasPath);
atlas = spm_read_vols(atlasStruct);

% Notify the user
fprintf('DONE: ATLAS data loaded.\n');



end
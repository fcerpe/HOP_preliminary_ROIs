function [refSpace, refStruct] = loadReference(opt, m, outputFolder)

% HOPLAB function to load a reference image used in resampling ROIs
%
% - takes method information (which refrence image to load)
% - if necessary, gunzips the image
% - copies image to output folder, renaming it
% - loads image and returns a struct element and a mars_space


% Extract parts of the filename
[~, refFilename, refExtension] = fileparts(m.referencePath);

% Create a new name for the reference image
newRefName = ['reference_task-', opt.task ,'_space-', opt.space, '.nii'];

% If needed, gunzip NIfTI file 
if strcmp(refExtension, '.gz')

    % Gunzip and overwrite referencePath with new one
    refPath = gunzip(m.referencePath, outputFolder);
    m.referencePath = refPath{1};

    % Get the filename without the extension
    atlasFilenameSplit = split(refFilename, '.');
    refFilename = atlasFilenameSplit{1};

elseif strcmp(refExtension, '.nii')

    % Copy the file from the source folder to the destination folder, and
    % re-name it (beta_0001.nii is not informative)
    % - forces overwrite
    newReferencePath = fullfile(outputFolder, newRefName);
    copyfile(m.referencePath, newReferencePath, 'f');
    m.referencePath = newReferencePath;

else
    error('No reference file found!');
end

% Load reference image for resampling

% Notify the user
fprintf('STEP: Getting REF data \n');

% Load the NIfTI file with SPM functions, use them to create marsbar_space
refSPM = spm_vol(m.referencePath);
refStruct = refSPM(1,1);
refSpace = mars_space(refSPM(1,1));

% Notify the user
fprintf('DONE: REF data loaded \n');


end
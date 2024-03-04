function [refSpace, refImgStruct] = loadReferenceImage(m, outputFolder)


% unzipping can be subfunction
% .nii if can be subfunction

[refPath, refFilename, refExt] = fileparts(m.referencePath);

if strcmp(refExt, '.gz')
    % Gunzip and extract gunzipped file path
    refPathCell = gunzipNiftiFile(m.referencePath, outputFolder);
    m.referencePath = refPathCell{1};

    % Get the filename without the .nii extension
    atlasFilenameSplit = split(refFilename, '.');
    refFilename = atlasFilenameSplit{1};

elseif strcmp(refExt, '.nii')
    % Copy the file from the source folder to the destination folder
    newrefImagePath = fullfile(outputFolder, strcat(refFilename, '.nii'));
    copyfile(m.referencePath, newrefImagePath);
    m.referencePath = newrefImagePath;
else
    error('No file ATLAS file found!');
end

% Load reference image for resampling
fprintf('STEP: Getting REF data \n');
refImgsStructLong = spm_vol(m.referencePath);
refImgStruct = refImgsStructLong(1,1);
refSpace = mars_space(refImgStruct);
fprintf('DONE: REF data loaded \n');


end
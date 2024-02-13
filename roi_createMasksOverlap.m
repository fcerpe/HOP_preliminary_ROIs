function [intersectedMask, outputFile] = roi_createMasksOverlap(specification, volumeDefiningImage, outputDir, saveImg)
% Adapted from bidspm/createROI
%
% See the ``demos/roi`` to see examples on how to use it.
%
% (C) Copyright 2021 CPP ROI developers

% Creates the first mask from neurosynth mask, then create the second
% ideally we want to loop over the masks and figure out
% if they are binary images or spheres...
nsImage = specification.mask1;
locaImage = specification.mask2;

isBinaryMask(nsImage);
isBinaryMask(locaImage);

nsMask = createRoi('mask', nsImage);
locaMask = createRoi('mask', locaImage);

locationsToSample = nsMask.global.XYZmm;

[~, nsMask.roi.XYZmm] = spm_ROI(locaMask, locationsToSample);

intersectedMask = setRoiSizeAndType(nsMask, 'intersection');

outputFile = [];

% Check that such intersection exists (i.e. there are voxels in the mask)
% otherwise skip the saving and throw a warning
if ~isempty(intersectedMask.roi.XYZmm)
    if saveImg
        outputFile = saveRoi(intersectedMask, volumeDefiningImage, outputDir);
    end

else
    warning('Intersection is empty. Will skip it')

end

end

function mask = setRoiSizeAndType(mask, type)
  mask.label = type;
  mask.descrip = type;
  mask.roi.size = size(mask.roi.XYZmm, 2);
end

function outputFile = saveRoi(mask, volumeDefiningImage, outputDir)

  hdr = spm_vol(volumeDefiningImage);
  if numel(hdr) > 1
    err.identifier =  'createRoi:not3DImage';
    err.message = sprintf(['the volumeDefininigImage:', ...
                           '\n\t%s\n', ...
                           'must be a 3D image. It seems to be 4D image with %i volume.'], ...
                          image, numel(hdr));
    error(err);
  end

  if ~strcmp(mask.def, 'sphere') && ...
      exist(mask.spec, 'file') == 2 && ...
      strcmp(spm_file(mask.spec, 'ext'), 'nii')
    checkRoiOrientation(volumeDefiningImage, mask.spec);
  end

  roiName = createRoiName(mask, volumeDefiningImage);

  % use the marsbar toolbox
  roiObject = maroi_pointlist(struct('XYZ', mask.roi.XYZmm, ...
                                     'mat', spm_get_space(volumeDefiningImage), ...
                                     'label', mask.label, ...
                                     'descrip', mask.descrip));

  % use Marsbar to save as a .mat and then convert that to an image
  % in the correct space
  outputFile = fullfile(outputDir, roiName);
  save_as_image(roiObject, outputFile);

  json = bids.derivatives_json(outputFile);
  bids.util.jsonencode(json.filename, json.content);

end

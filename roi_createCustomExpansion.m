function [mask, outputFile, report] = roi_createCustomExpansion(specification, VolumeDefiningImage, OutputDir, saveImg, attempt)
%
% Returns a mask to be used as a ROI by ``spm_summarize``.
% Can also save the ROI as binary image.
% See 'createRoi for more information.
%
% Customization by Filippo Cerpelloni,
% orginal script (C) Copyright 2023 CPP ROI developers
%
% Customization consists of two modifications. This script additionally
% checks
%
% 1. that the expansion of the sphere actually corresponds to an increase
% in the voxel size. If not, it's an indication that the cluster from which
% we started drawing the sphere has been mapped fully and there are no more
% voxels
%
% 2. that the radius is smaller than 15mm. If not, we can safely assume
% that we are drawing an ROI that is too large, given the nature of our
% original ROIs
%
% If one of this conditions are not met, we iterate again on a more lax
% threshold

%% Compute ROI

% First iteration parameters: 
% - the initial type of ROI is an expansion
type = 'expand';

% - we are not in a subsequent attempt
redo = false;


switch type

    case 'sphere'

        sphere = specification;

        mask.def = type;
        mask.spec = sphere.radius;
        mask.xyz = sphere.location;

        if size(mask.xyz, 1) ~= 3
            mask.xyz = mask.xyz';
        end

        mask = spm_ROI(mask);
        mask.roi.XYZmm = [];

        mask = createRoiLabel(mask);

    case 'mask'

        roiImage = specification;

        isBinaryMask(roiImage);

        mask = struct('XYZmm', []);
        mask = defineGlobalSearchSpace(mask, roiImage);

        % in real world coordinates
        mask.global.XYZmm = returnXYZm(mask.global.hdr.mat, mask.global.XYZ);

        assert(size(mask.global.XYZmm, 2) == sum(mask.global.img(:)));

        locationsToSample = mask.global.XYZmm;

        mask.def = type;
        mask.spec = roiImage;
        [~, mask.roi.XYZmm, j] = spm_ROI(mask, locationsToSample);

        mask.roi.XYZ = mask.global.XYZ(:, j);

        mask = setRoiSizeAndType(mask, type);

        mask = createRoiLabel(mask);

    case 'intersection'

        % load the thresholded mask corresponding to the iteration nb
        switch attempt
            case 1, roiImage = specification.mask1; % p = 0.001
            case 2, roiImage = specification.mask2; % p = 0.01
            case 3, roiImage = specification.mask3; % p = 0.05
            case 4, roiImage = specification.mask4; % p = 0.1
        end

        % load and check mask
        sphere = specification.maskSphere;

        isBinaryMask(roiImage);

        mask = createRoi('mask', roiImage);
        mask2 = createRoi('sphere', sphere);

        locationsToSample = mask.global.XYZmm;

        [~, mask.roi.XYZmm] = spm_ROI(mask2, locationsToSample);

        mask = setRoiSizeAndType(mask, type);

        mask = createRoiLabel(mask);

    case 'expand'

        % load the thresholded mask corresponding to the iteration nb
        switch attempt
            case 1, roiImage = specification.mask1; % p = 0.001
            case 2, roiImage = specification.mask2; % p = 0.01
            case 3, roiImage = specification.mask3; % p = 0.05
            case 4, roiImage = specification.mask4; % p = 0.1
        end

        % If all attemptes go wrong, the problem is likely before the ROI
        % creation. Notify the user
        if attempt == 5
            error(['Even with progressively lax thresholds, there are not enough voxels meeting our criteria.',...
                   'Check the previous steps of the pipeline for possible errors.']);
        end

        sphere = specification.maskSphere;

        isBinaryMask(roiImage);

        % check that input image has at least enough voxels to include
        maskVol = spm_read_vols(spm_vol(roiImage));
        totalNbVoxels = sum(maskVol(:));

        if sphere.maxNbVoxels > totalNbVoxels
            % notify the user
            fprintf('\nNumber of voxels requested greater than the total number of voxels in this mask\n');

            % Go to the next iteration and keep track of the nb of attempts
            redo = true;
            attempt = attempt+1;

        else
            spec  = struct('mask1', roiImage, ...
                           'mask2', sphere);
    
            % take as radius step the smallest voxel dimension of the roi image
            hdr = spm_vol(roiImage);
            dim = diag(hdr.mat);
            radiusStep = min(abs(dim(1:3)));
    
            % determine maximum radius to expand to
            maxRadius = hdr.dim .* dim(1:3)';
            maxRadius = max(abs(maxRadius));
    
            fprintf(1, '\n Expansion:');
    
            previousSize = 0;
    
            while  true
                mask = createRoi('intersection', spec);
                mask.roi.radius = spec.mask2.radius;
    
                fprintf(1, '\n radius: %0.2f mm; roi size: %i voxels', ...
                    mask.roi.radius, ...
                    mask.roi.size);
    
                % Check that size is indeed expanding, if not go to the
                % next attempt
                if mask.roi.size == previousSize
                    fprintf(['\n An increase in the size did not increase the number of voxels. ' ...
                             'This cluster''s size is not sufficient.\n' ...
                             'Will use lower threshold\n']);
                    redo = true;
                    attempt = attempt+1;
                    break
                end
    
                % Check size of the sphere, if too large go to the next
                % attempt
                if mask.roi.radius >= 15
                    fprintf(['\n Radius too large\nIt would not make sense for our area ' ...
                             'to be this large and still not have enough voxels.\n' ...
                             'Will use lower threshold\n']);
                    redo = true;
                    attempt = attempt+1;
                    break
                end    
    
                % if there are no more voxels in the whole mask, stop
                if mask.roi.size > sphere.maxNbVoxels
                    break
                end
    
                % Similarly, stop if the radius is larger than the mask itself
                if mask.roi.radius > maxRadius
                    error('sphere expanded beyond the dimension of the mask.');
                end
    
                spec.mask2.radius = spec.mask2.radius + radiusStep;
                previousSize = mask.roi.size;
            end
        end
        
        if redo
            % At any point, the current expansion process was not
            % sufficient to extract an ROI. 
            % Re-iterate with a different threshold
            %
            % (this scripting is dangerous, because it calls itself. use carefully)
            mask = roi_createCustomExpansion(specification, VolumeDefiningImage, OutputDir, saveImg, attempt);
        end


        fprintf(1, '\n');

        % create label for the ROI
        mask.xyz = sphere.location;
        mask = setRoiSizeAndType(mask, type);
        mask = createRoiLabel(mask);
end

%% save ROI

outputFile = [];

if saveImg
    outputFile = saveRoi(mask, VolumeDefiningImage, OutputDir);
end

% additionally save a report of the extraction process: 
% which threshold was used, how many voxels are in the final mask, how
% large was the final sphere radius
if not(redo)

    % Directly write the document

    % open report - with date notation to have a new one each time we try
    reportFile = dir(['reports/roi_reports_expansion_' date '.txt']);
    if not(isempty(reportFile))
        report = readcell(reportFile.name);
    else
        report = [];
    end

    % add a new line
    % - method
    reportMethod = mask.def;

    % - area
    splitLoca = split(specification.mask1,{'desc-','_label'});
    if startsWith(splitLoca{2},'french'), reportArea = 'VWFA';
    else, reportArea = 'LO';
    end

    % - mask threshold
    switch attempt
        case 1, reportPvalue = '0.001';
        case 2, reportPvalue = '0.01';
        case 3, reportPvalue = '0.05';
        case 4, reportPvalue = '0.1';
    end

    % - attempt nb
    reportAttempt = attempt;

    % - sphere radius
    reportRadius = mask.roi.radius;

    % - nb of voxels
    reportVoxels = mask.roi.size;

    % Add everything to the new row
    report = vertcat(report, {reportArea, reportMethod, reportAttempt, reportPvalue, reportRadius, reportVoxels});

    % save report
    writecell(report,['reports/roi_report_expansions_' date '.txt'])
end

end


%% Support functions

function mask = defineGlobalSearchSpace(mask, image)
% gets the X, Y, Z subscripts of the voxels included in the ROI

mask.global.hdr = spm_vol(image);
mask.global.img = logical(spm_read_vols(mask.global.hdr));

[X, Y, Z] = ind2sub(size(mask.global.img), find(mask.global.img));

% XYZ format
mask.global.XYZ = [X'; Y'; Z'];
mask.global.size = size(mask.global.XYZ, 2);

end

function XYZmm = returnXYZm(transformationMatrix, XYZ)
% apply voxel to world transformation

XYZmm = transformationMatrix(1:3, :) * [XYZ; ones(1, size(XYZ, 2))];
end

function mask = setRoiSizeAndType(mask, type)
mask.def = type;
mask.roi.size = size(mask.roi.XYZmm, 2);
end

%% Re-adapted function to create labels
function mask = createRoiLabel(mask)

switch mask.def

    case 'sphere'

        mask.label = sprintf('sphere%0.0fx%0.0fy%0.0fz%0.0f', ...
            mask.spec, ...
            mask.xyz);

        % change any minus coordinate (x = -67) to minus (xMinus67)
        % SUPER ugly but any minus will mess up the bids parsing otherwise
        mask.label = strrep(mask.label, '-', 'Minus');

        mask.descrip = sprintf('%s at [%0.1f %0.1f %0.1f]', ...
            mask.str, ...
            mask.xyz);

    case 'expand'

        mask.label = sprintf('%sVox%i', ...
            mask.def, ...
            mask.roi.size);

        mask.descrip = sprintf('%s from [%0.0f %0.0f %0.0f] till %i voxels', ...
            mask.def, ...
            mask.xyz, ...
            mask.roi.size);

    otherwise

        mask.label = mask.def;
        mask.descrip = mask.def;

end

end

%% Re-adapted funcion to save the rois
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

if strcmp(mask.def, 'sphere')

    [~, mask.roi.XYZmm] = spm_ROI(mask, volumeDefiningImage);
    mask = setRoiSizeAndType(mask, mask.def);

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

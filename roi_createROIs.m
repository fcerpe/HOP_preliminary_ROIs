%% Create ROIs based on peak coordinates from subjects' localizers
%
% From individual activation peaks identified by the localizers, extract
% the expansion for a specified number of voxels, then overlap the
% resulting mask with those extracted from neurosynth.org and relative to
% the same contrast: in visbra, 'visual words' for VWFA and 'objects' for LO
% 
% Outputs:
% - ROIs (duh) in the folder specified in roi_option
% - expansionReport_date.txt with a recap of the expansion process. It will
%   show 'area', 'method', 'iteration', 'p', 'sphere radius', 'number of
%   voxels'
%
% Steps:
% 1) get peak coordinates - requires hard-coded peaks, see roi_getMNIcoords
%
% 2) draw expansion with different thresholds 
%    - requires already-made contrasts with p = 0.001; 0.01; 0.05; 0.1
%      (made with stats_localizer)
% 
% 3) overlap localizer expansion with neurosynth thresholded mask
% 
% 4) resolve potential overlaps between areas (momentarily in a second
%    script)

%% Get the ROIs 

% get the resliced masks in the folder. Only 'visualWords' and 'objects'
% fit the criterium
neurosynthMasks = dir("masks/neurosynth/r*");

% Get the peaks for each sub / area
mni = roi_getMNIcoords(opt.subjects);

% Get information over the ROIs
% - names: VWFA, lLO, rLO
roiNames = opt.roiList;

% - contrasts to use:
contrastName = {'frenchGtScrambled', 'brailleGtScrambled', 'drawingGtScrambled', 'drawingGtScrambled'};

% - hemisphere: 
hemiName = {'L','L','L','R'};

% - folder contating the localizer activations:
localizerStatsFolder = ['task-visualLocalizer_space-',opt.space{1},'_FWHM-6_node-localizerGLM'];


% specify the number of voxels to keep
for iSub = 1:length(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];

    % look for folder in cpp_spm-rois
    % if does not exists, create it
    if isempty(dir(fullfile(opt.dir.rois, subName)))
        mkdir(fullfile(opt.dir.rois, subName))
    end

    % for each region this subject has
    for iReg = 1:length(mni{1}(:,1))

        % Skip undefinded regions and braille contrast in VWFA
        if all(not(isnan(mni{iSub}(iReg, :)))) && iReg ~= 2 

            %% STEP 1 : expansion around peak coordinates
            % Expansion takes the localizer activation, a peak coordinate
            % and expands form the coordinate creating progressively large
            % spheres. Those spheres are then intersected with the
            % localizer mask to extrat only the relevant voxels for a given
            % contrast

            % Take the peak of activation
            roiCenter = mni{iSub}(iReg, :);
            
            % Get a reference image for reslicing
            dataImage = fullfile(opt.dir.stats, subName, localizerStatsFolder, 'beta_0001.nii');

            % Get the filename of the corresponding contrasts with
            % different thresholds
            mask001InDir = dir(fullfile(opt.dir.stats, subName, localizerStatsFolder, ...
                [subName, '_task-visualLocalizer_space-',opt.space{1},'_desc-', contrastName{iReg} ,'*_p-0pt001_k-0_MC-none_mask.nii']));
            mask01InDir = dir(fullfile(opt.dir.stats, subName, localizerStatsFolder, ...
                [subName, '_task-visualLocalizer_space-',opt.space{1},'_desc-', contrastName{iReg} ,'*_p-0pt010_k-0_MC-none_mask.nii']));
            mask05InDir = dir(fullfile(opt.dir.stats, subName, localizerStatsFolder, ...
                [subName, '_task-visualLocalizer_space-',opt.space{1},'_desc-', contrastName{iReg} ,'*_p-0pt050_k-0_MC-none_mask.nii']));
            mask1InDir = dir(fullfile(opt.dir.stats, subName, localizerStatsFolder, ...
                [subName, '_task-visualLocalizer_space-',opt.space{1},'_desc-', contrastName{iReg} ,'*_p-0pt100_k-0_MC-none_mask.nii']));
            
            % Get the full pathof each thresholded contrast (if exists)
            localizer001Mask = fullfile(mask001InDir.folder, mask001InDir.name);

            if not(isempty(mask05InDir))
                localizer05Mask = fullfile(mask05InDir.folder, mask05InDir.name);
            end

            if not(isempty(mask01InDir))
                localizer01Mask = fullfile(mask01InDir.folder, mask01InDir.name);
            end

            if not(isempty(mask1InDir))
                localizer1Mask = fullfile(mask1InDir.folder, mask1InDir.name);
            end

            % specify the sphere parameters for each of them
            sphereParams = struct;
            sphereParams.location = roiCenter;
            sphereParams.radius = 1; % starting radius
            sphereParams.maxNbVoxels = opt.numVoxels;

            % Add all the masks that are present in the localizer folder
            specification = struct('mask1', localizer001Mask, 'maskSphere', sphereParams);
            if not(isempty(mask01InDir))
                specification.mask2 = localizer01Mask;
            end
            if not(isempty(mask05InDir))
                specification.mask3 = localizer05Mask;
            end
            if not(isempty(mask1InDir))
                specification.mask4 = localizer1Mask;
            end

            % specify the path for each subject
            outputPath = fullfile(opt.dir.rois, subName);

            % Notify the user
            fprintf('\nWorking on %s, area %s \n',subName, char(roiNames{iReg}));
           
            % Compute expansion 
            [~, sphereMaskName] = roi_createCustomExpansion(specification, dataImage, outputPath, opt.saveROI, 1);

            % modify the name:
            % from mentions of the localizer mask to a more bids-like name
            sphereMaskNameOnly = sphereMaskName(1:end-4);
            findNbVox = split(sphereMaskNameOnly, {'Vox','_desc'});

            % Custom name
            bidslikeName = fullfile(opt.dir.rois, subName, [subName,'_hemi-',hemiName{iReg},'_space-',opt.space{1}, ...
                                                            '_label-',char(roiNames{iReg}),'_voxels-',findNbVox{2},'_mask']);

            % Rename .nii and .json files
            movefile(sphereMaskName, [bidslikeName,'.nii'],'f')
            movefile([sphereMaskNameOnly,'.json'], [bidslikeName,'.json'],'f')

            % Reslice the mask to match our voxel size and space
            sphereMask = resliceRoiImages(dataImage, [bidslikeName, '.nii']);


            %% STEP 2 : intersection with neurosynth mask
            % Overlap the newly-created mask to the neurosynth mask to ensure 
            % accurate localization.
            % Assumption is that expansion method grows indiscriminately
            % around peak activation. Overlapping with a neurosynth mask
            % should make sure that the ROI is bounded by neuroscientific
            % constrains. 
            % (could be more efficient)

            % Get the neurosynth masks
            % - 'objects' is first because of [B]ilateral hemisphere 
            % - 'visualWords' has only [L]eft hemisphere
            switch iReg
                case {1,2}, nsMask = fullfile(neurosynthMasks(2).folder, neurosynthMasks(2).name);
                case {3,4}, nsMask = fullfile(neurosynthMasks(1).folder, neurosynthMasks(1).name);
            end

            % Load and cast the nifti file as 'uint8' to be read properly
            temp = load_nii(nsMask);
            temp.img = cast(temp.img, 'uint8');
            save_nii(temp, nsMask);

            % specify the objects to pass: mask + sphere or mask + mask
            specMasks = struct('mask1', nsMask, ...
                               'mask2', sphereMask);

            % Compute intersection between expansion and neurosynth mask
            [intersectedMask, intersectedName] = roi_createMasksOverlap(specMasks, dataImage, outputPath, opt.saveROI);

            % Remove file extension from name to be used in renaming
            intersectedNameOnly = intersectedName(1:end-4);

            % Create custom new name
            intersectedNewName = fullfile(opt.dir.rois, subName, [subName, '_hemi-' hemiName{iReg}, ...
                                                                  '_space-', opt.space{1}, ...
                                                                  '_atlas-neurosynth_method-expansionIntersection_label-', ...
                                                                  roiNames{iReg}, '_mask']);
 
            % Rename intersected mask
            movefile(intersectedName, [intersectedNewName,'.nii'],'f')
            movefile([intersectedNameOnly,'.json'], [intersectedNewName,'.json'],'f')

            % reslice the masks
            intersectedMask = resliceRoiImages(dataImage, [intersectedNewName, '.nii']);

        end
    end
end







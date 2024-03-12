function hop_roi_expansion(opt, m)

% FUNCTION TO GENERATE EXPANSION ROIS IN NIFTI FORMAT
%
% This script generates ROIs (Regions of Interest) starting from peak coordinates
% and a activation map, defined in 'hop_option'. 
% It creates progressively larger sphere ROIs around the peak coordinates 
% provided, until a certain target number of voxels if reached. If the
% sphere cannot expand further, it will relax the threshold until p < 0.05
% uncorr. The resulting ROIs are saved as NIfTI files in an out folder.
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
%     - iterates over each speficied radius
%     - if not present, creates an output folder for each radius
%     - iterates over each ROI requested
%     - creates sphere ROIs for the left and right hemispheres using the specified MNI coordinates
%     - resamples the ROIs to match the reference image space
%     - verifies the compatibility of the resampled ROI and the reference image
%     - stores the resampled ROIs in a cell array
%     - merges the resampled ROIs if coordinates are provided for both hemispheres
%     - saves the merged ROI as a NIfTI file in the output folder
%     - copies the script file to the output folder for replicability.
%
% The script is executed for each radius specified in the radii vector%
% Dependencies: MarsBaR, SPM
%
% Author: Filippo Cerpelloni
% Date: March 2024

% iterate over every subject
for iSub = 1:length(m.subjects)

    % Get subject name
    subName = ['sub-', m.subjects{iSub}];
    fprintf('\n###### Processing %s ######\n', subName);

     % Iterate for each area requested
    for iROI = 1:length(m.roisToCreate)

        % Take the current ROI that is being extracted
        currROI = m.roisToCreate(iROI);

        % Create output folder
        % skip method details, in this case uses the target nb of voxels
        % already present in 'm'
        outputFolder = createOutputFolder(opt, m, [], subName);

        % Notify the user
        fprintf('\n--- Processing ROI: %s -  %s ---\n', currROI.area, currROI.contrast);

        % Get the subject's folder
        % Does not check for presence of the task, node, etc.
        % Refer to 'checkSubjectHasTask' function to improve
        subjectFolder = fullfile(opt.dir.stats, subName, ...
                                ['task-', opt.task,'_space-',opt.space,'_FWHM-', num2str(opt.fwhm),'_node-',opt.node]); 

        % Load SPM model 
        [model, loadSuccess] = loadSPMmodel(subjectFolder, subName, currROI.task); % Import the SPM model file

        % Skip this contrast / ROI if the SPM model fails to load
        if ~loadSuccess, continue 
        end
        
        % Load T constrast and apply threshold (voxel with p < .001)
        [tConImg, tStruct] = loadTContrastImage(model, currROI.contrast, subjectFolder);
        [thresholdedImage, finalPThresh] = calculateAndApplyThreshold(tConImg, model, [], [], m.targetNbVoxels);




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











end
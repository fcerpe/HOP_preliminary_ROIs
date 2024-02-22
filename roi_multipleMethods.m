%% Main script to create ROIs in different ways, for comparison
%
% To settle a debate between Hans and Olivier: which ROI definition method
% should we use?
% 1) 10mm sphere around neurosynth coordinates
%    VWFA = -44 -56 -16, keyword: visual words
%    lLO  = -47 -70 -5
%    rLO  =  47 -70 -5, keyword: objects
%
% 2) 10mm sphere around individual coordinates (extracted from a 10mm
%    sphere around the neurosynth coordinates)
%
% 3) 8mm sphere around individual coordinates
%
% 4) expanding volume of 150 voxels around the individual coordinates
%    not done in this point. For now we apply the binary mask to it
%
% 5) anatomically-defined ROIs


%% Create / initialize the basics
clear;
clc;

% add bidspm and init it
addpath '/lib/bidspm'
bidspm;

% get options
opt = roi_option();

%% Method #1 - Take everything that is where the area is supposed to be

% get general coordinates (manual input)
mni{1}(1, 1:3) = [-4.420000e+01, -5.480000e+01, -1.540000e+01];  % VWFA FRE
mni{1}(3, 1:3) = [-4.680000e+01, -7.040000e+01, -0.500000e+01];  %  LO LEFT
mni{1}(5, 1:3) = [ 4.680000e+01, -7.040000e+01, -0.500000e+01];  %  LO RIGHT

roiNames = opt.roiList;
save('ROI-trial_neurosynth-coords.mat', 'roiNames', 'mni');

% Get the ROIs (actually just the spheres)
for iSub = 1:length(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];

    % for each region this subject has
    for iReg = 1:length(mni{1}(:,1))

        % if the region is defined (vwfa-br for all, but in some is not present)
        if not(isnan(mni{1}(iReg, :)))

            % Get the center
            ROI_center = mni{1}(iReg, :);

            % Get the name of the roi for filename
            switch iReg
                case 1, regName = 'VWFAfr';
                case 2, regName = 'VWFAbr';
                case 3, regName = 'lLO';
                case 4, regName = 'lpFS';
                case 5, regName = 'rLO';
                case 6, regName = 'rpFS';
            end

            % Set up bids-like name
            bidslikeName = [subName,'_','space-MNI','_', 'trial-NeurosynthCoords_' , ...
                'label-',regName,'_','radius-',num2str(opt.radius),'mm','_mask'];

            betaReference = fullfile(opt.dir.stats, subName, ...
                'task-wordsDecoding_space-IXI549Space_FWHM-2', ...
                'beta_0001.nii');

            % specify the sphere characteristics for each of them
            sphereParams = struct;
            sphereParams.location = ROI_center;
            sphereParams.radius = 10;

            % specify the path for each subject
            outputPath = [opt.dir.rois,'/',subName];

            sphereMask = createRoi('sphere', sphereParams, betaReference, outputPath, opt.saveROI);

            path = '../../outputs/derivatives/cpp_spm-rois/';
            nativeName = ['label-',sphereMask.label,'_mask'];
            movefile([path, subName, '/', nativeName,'.nii'], [path, subName, '/', bidslikeName,'.nii'],'f')
            movefile([path, subName, '/', nativeName,'.json'], [path, subName, '/', bidslikeName,'.json'],'f')

            % reslice
            roiPath = fullfile(opt.dir.rois, subName, [bidslikeName, '.nii']);
            sphereMask = resliceRoiImages(betaReference,roiPath);

        end
    end
end


%% Method #2 and #3 - 8mm and 10mm spheres around individual coordinates (one script for efficiency)
% Technically intersections between sphere and brain mask

% get individual coordinates
mni = roi_getMNIcoords(opt.subjects);

roiNames = opt.roiList;
save('ROI-trial_indivudal-coords.mat', 'roiNames', 'mni');

% Get the ROIs (actually just the spheres)
for iSub = 1:length(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];

    % for each region this subject has
    for iReg = 1:length(mni{1}(:,1))

        % if the region is defined (vwfa-br for all, but in some is not present)
        if not(isnan(mni{iSub}(iReg, :)))
            % Get the center
            ROI_center = mni{iSub}(iReg, :);

            % Get the name of the roi for filename
            switch iReg
                case 1, regName = 'VWFAfr';
                case 2, regName = 'VWFAbr';
                case 3, regName = 'lLO';
                case 4, regName = 'lpFS';
                case 5, regName = 'rLO';
                case 6, regName = 'rpFS';
            end

            for rad = 8:2:10
                % Set up bids-like name
                bidslikeName = [subName,'_','space-MNI','_', 'trial-IndividualCoords_' , ...
                    'label-',regName,'_','radius-',num2str(rad),'mm','_mask'];

                betaReference = fullfile(opt.dir.stats, subName, 'task-wordsDecoding_space-IXI549Space_FWHM-2', 'beta_0001.nii');

                % anatomical mask is brain mask, make sure that the sphere
                % doesn't excedd the limits of what makes sense
                brainMask = fullfile(opt.dir.stats, subName, 'task-wordsDecoding_space-IXI549Space_FWHM-2', 'mask.nii'); 

                % specify the sphere characteristics for each of them
                sphereParams = struct;
                sphereParams.location = ROI_center;
                sphereParams.radius = rad;

                % % specify the object to pass: mask + sphere
                specification = struct('mask1', brainMask, ...
                                       'mask2', sphereParams);

                % specify the path for each subject
                outputPath = [opt.dir.rois,'/',subName];

                sphereMask = createRoi('intersection', specification, betaReference, outputPath, opt.saveROI);

                path = '../../outputs/derivatives/cpp_spm-rois/';
                nativeName = ['label-',sphereMask.label,'_mask'];
                movefile([path, subName, '/', nativeName,'.nii'], [path, subName, '/', bidslikeName,'.nii'],'f')
                movefile([path, subName, '/', nativeName,'.json'], [path, subName, '/', bidslikeName,'.json'],'f')

                % reslice
                roiPath = fullfile(opt.dir.rois, subName, [bidslikeName, '.nii']);
                sphereMask = resliceRoiImages(betaReference,roiPath);
            end
        end
    end
end

%% Method #4 - expanding volume of 100vx starting from individual peaks

% start with a new report
delete expansionReport.txt

% get the resliced masks in the folder. Only 'visualWords' and 'objects'
% fit the criterium
neurosynthMasks = dir("masks/r*");

mni = roi_getMNIcoords(opt.subjects);

roiNames = opt.roiList;
save('ROI_trial-expansion.mat', 'roiNames', 'mni');

% which name to assign: 1 to 6 based on the region
regName = {'VWFAfr','VWFAbr','lLO','rLO'};
contrastName = {'frenchGtScrambled', 'brailleGtScrambled', 'drawingGtScrambled', 'drawingGtScrambled'};
hemiName = {'L','L','L','R'};

% specify th number of voxels to keep
opt.numVoxels = 108;

for iSub = 1:length(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];

    % for each region this subject has
    for iReg = 1:length(mni{1}(:,1))

        % if the region is defined (a subject does not have a
        % braille-defined VWFA at all)
        % AND if the region is not the braille-defined VWFA (not useful for
        % analysis at this stage)
        if all(not(isnan(mni{iSub}(iReg, :)))) && iReg ~= 2 

            % STEP 1 : expansion around peak coordinates

            % Get the center
            ROI_center = mni{iSub}(iReg, :);
            
            % Get the reference image
            dataImage = fullfile(opt.dir.stats, subName, 'task-wordsDecoding_space-IXI549Space_FWHM-2', 'beta_0001.nii');

            % Get the filename of the corresponding contrast
            mask001InDir = dir(fullfile(opt.dir.stats, subName, 'task-visualLocalizer_space-IXI549Space_FWHM-6', ...
                [subName, '_task-visualLocalizer_space-IXI549Space_desc-', contrastName{iReg} ,'_*_p-0pt001_k-0_MC-none_mask.nii']));
            mask01InDir = dir(fullfile(opt.dir.stats, subName, 'task-visualLocalizer_space-IXI549Space_FWHM-6', ...
                [subName, '_task-visualLocalizer_space-IXI549Space_desc-', contrastName{iReg} ,'_*_p-0pt010_k-0_MC-none_mask.nii']));
            mask05InDir = dir(fullfile(opt.dir.stats, subName, 'task-visualLocalizer_space-IXI549Space_FWHM-6', ...
                [subName, '_task-visualLocalizer_space-IXI549Space_desc-', contrastName{iReg} ,'_*_p-0pt050_k-0_MC-none_mask.nii']));
            mask1InDir = dir(fullfile(opt.dir.stats, subName, 'task-visualLocalizer_space-IXI549Space_FWHM-6', ...
                [subName, '_task-visualLocalizer_space-IXI549Space_desc-', contrastName{iReg} ,'_*_p-0pt100_k-0_MC-none_mask.nii']));
            
            % Get the full path: folder + name
            localizer001Mask = fullfile(mask001InDir.folder, mask001InDir.name);
            if not(isempty(mask05InDir)), localizer05Mask = fullfile(mask05InDir.folder, mask05InDir.name);
            end
            if not(isempty(mask01InDir)), localizer01Mask = fullfile(mask01InDir.folder, mask01InDir.name);
            end
            if not(isempty(mask1InDir)),  localizer1Mask = fullfile(mask1InDir.folder, mask1InDir.name);
            end

            % specify the sphere characteristics for each of them
            sphereParams = struct;
            sphereParams.location = ROI_center;
            sphereParams.radius = 1; % starting radius
            sphereParams.maxNbVoxels = opt.numVoxels;

            specification = struct('mask1', localizer001Mask, 'maskSphere', sphereParams);

            % add other masks, if not empty
            if not(isempty(mask01InDir)),   specification.mask2 = localizer01Mask;
            end
            if not(isempty(mask05InDir)),   specification.mask3 = localizer05Mask;
            end
            if not(isempty(mask1InDir)),    specification.mask4 = localizer1Mask;
            else

            end
            % specify the path for each subject
            outputPath = [opt.dir.rois,'/',subName];

            fprintf('\nWorking on %s, area %s \n',subName, char(regName{iReg}));
           
            % expansion based on localizer mask(s): 
            % - draw progressively larger spheres around the peak
            %   coordinates and take only the active voxels in the mask
            % - if the sphere expands without including other voxels, means
            %   the cluster is finished: use a more lax threshold and try
            %   again until we reach our number of voxels
            [sphereMask, sphMaskName] = roi_createCustomExpansion(specification, dataImage, outputPath, opt.saveROI, 1);

            % manipulations on the filename to make it more readable
            sphMaskJustTheName = sphMaskName(1:end-4);
            findVox = split(sphMaskJustTheName, {'Vox','_desc'});

            % Set up a custom bids-like name
            bidslikeName = fullfile(opt.dir.rois, subName, [subName,'_','space-MNI','_', 'trial-Expansion_' , ...
                'label-',char(regName{iReg}),'_','voxels-',findVox{2},'vx','_mask']);

            % rename both the .nii and the .json files to the new name
            movefile(sphMaskName, [bidslikeName,'.nii'],'f')
            movefile([sphMaskJustTheName,'.json'], [bidslikeName,'.json'],'f')

            % reslice the mask to make sure fits our voxel size and voxel
            % space
            sphereMask = resliceRoiImages(dataImage, [bidslikeName, '.nii']);


            % STEP 2 : intersection with neurosynth mask
            % (could be written to be more efficient, at the moment
            % (20/02/23) it's a WIP)

            % Get the neurosynth masks
            % 'objects' is first because of [B]ilateral hemisphere 
            % 'visualWords' has only [L]eft hemisphere
            switch iReg
                case {1,2}, nsMask = fullfile(neurosynthMasks(2).folder, neurosynthMasks(2).name);
                case {3,4}, nsMask = fullfile(neurosynthMasks(1).folder, neurosynthMasks(1).name);
            end

            % Open the image to make sure it's uint8
            % otherwise it's not read properly
            temp = load_nii(nsMask);
            temp.img = cast(temp.img, 'uint8');
            save_nii(temp, nsMask);

            % specify the objects to pass: mask + sphere or mask + mask
            specMasks = struct('mask1', nsMask, ...
                               'mask2', sphereMask);

            % specify the path for each subject
            outputPath = [opt.dir.rois,'/',subName];

            % intersection between expansion and neurosynth mask: 
            [intersectedMask, intersectedName] = roi_createMasksOverlap(specMasks, dataImage, outputPath, opt.saveROI);

            % Rename .json and .nii files of both masks to have more readable names
            % Remove file extension from name
            intersectedJustTheName = intersectedName(1:end-4);

            % New names
            intersectedNewName = fullfile(opt.dir.rois, subName, ...
                                [subName, '_hemi-' hemiName{iReg} ...
                                '_space-MNI_atlas-neurosynth_method-expansionIntersection_label-' regName{iReg} '_mask']);
 
            % Rename intersection
            movefile(intersectedName, [intersectedNewName,'.nii'],'f')
            movefile([intersectedJustTheName,'.json'], [intersectedNewName,'.json'],'f')

            % reslice the masks
            intersectedMask = resliceRoiImages(dataImage, [intersectedNewName, '.nii']);

        end
    end
end


%% #5 Anatomically-defined ROIs
% ROIs created using SPM Anatomy toolbox
% (https://github.com/inm7/jubrain-anatomy-toolbox)
%
% CPP_ROI has 'getDataFromIntersection' to extract an ROI from a mask and a
% sphere

% which name to assign: 1 to 6 based on the region
regName = {'VWFAfr','lLO','rLO'};

mni = roi_getMNIcoords(opt.subjects);

roiNames = opt.roiList;
save('ROI-trial_indivudal-coords.mat', 'roiNames', 'mni');

% Get the ROIs (actually just the spheres)
for iSub = 1:length(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];

    specificRegions = [1 2 4 6];

    % for each region this subject has
    for iReg = 1:length(mni{1}(:,1))

        % if the region is defined (vwfa-br for all, but in some is not present)
        % and if there is a definitve anatomical mask: applicable to 
        % VWFAs and lpFS
        if ismember(iReg,specificRegions) && not(all(isnan(mni{iSub}(iReg, :))))

            % Sphere data: get the center, radius is fixed
            ROI_center = mni{iSub}(iReg, :);
            opt.radius = 10;

            % Set up bids-like name to save ROI
            bidslikeName = [subName,'_','space-MNI','_','trial-AnatIntersection_', ...
                            'label-',regName{iReg},'_','radius-',num2str(opt.radius),'mm','_mask'];

            % Get the reference, a random beta
            betaReference = fullfile(opt.dir.stats, subName, 'task-wordsDecoding_space-IXI549Space_FWHM-2', 'beta_0001.nii');

            % Get the anatomical mask based on the contrast
            % left fg for vwfa
            if ismember(iReg,[1 2 4])
                anatomicalMask = 'JuBrain_ROIs/sub-none_space-MNI_desc-leftFG_label-0000_mask.nii'; 
            else % right hemisphere
                anatomicalMask = 'JuBrain_ROIs/sub-none_space-MNI_desc-rightFG_label-0000_mask.nii';
            end

            % Reslice it for compatibility with our 2.6mm voxel size
            anatomicalMask = resliceRoiImages(betaReference,anatomicalMask);

            % convert the data-type of the nii file from 'single' to
            % 'uint8', otherwise isBinaryMask does not work
            thisMask = load_nii(anatomicalMask);
            thisMask.img = cast(thisMask.img,'uint8');
            save_nii(thisMask,anatomicalMask);

            % specify the sphere characteristics for each of them
            sphereParams = struct;
            sphereParams.location = ROI_center;
            sphereParams.radius = opt.radius;

            % specify the object to pass: mask + sphere
            specification = struct('mask1', anatomicalMask, ...
                                   'mask2', sphereParams);

            % specify the path for each subject
            outputPath = [opt.dir.rois,'/',subName];

            sphereMask = createRoi('intersection', specification, betaReference, outputPath, opt.saveROI);

            path = '../../outputs/derivatives/cpp_spm-rois/';
            if iReg == 6
                nativeName = 'rsub-none_space-MNI_label-0000Intersection_desc-rightFG_mask';                
            else
                nativeName = 'rsub-none_space-MNI_label-0000Intersection_desc-leftFG_mask';
            end
            movefile([path, subName, '/', nativeName,'.nii'], [path, subName, '/', bidslikeName,'.nii'],'f')
            movefile([path, subName, '/', nativeName,'.json'], [path, subName, '/', bidslikeName,'.json'],'f')

            % reslice
            roiPath = fullfile(opt.dir.rois, subName, [bidslikeName, '.nii']);
            sphereMask = resliceRoiImages(betaReference,roiPath);

        end
    end
end

%% Method #6 - Neurosynth mask
% Mask from neurosynth: https://www.neurosynth.org/analyses/terms/words/
% Edited in bspmview to get a threshold of 7 (arbitrary value)

% masks to try quickly: diffrent arbirtary thresholds to see sub-007 and
% sub-008 (limit cases) 
neurosynthMasks = {'rhemi-L_space-MNI_atlas-neurosynth_label-visualWords_thresh-5_voxels-2557_mask.nii', ...
                   'rhemi-L_space-MNI_atlas-neurosynth_label-visualWords_thresh-6_voxels-1768_mask.nii', ...
                   'rhemi-L_space-MNI_atlas-neurosynth_label-visualWords_thresh-7_voxels-1222_mask.nii'};

mni = roi_getMNIcoords(opt.subjects);

roiNames = opt.roiList;
save('ROI-trial_indivudal-coords.mat', 'roiNames', 'mni');

% which name to assign: 1 to 6 based on the region
regName = {'VWFAfr'}; %,'VWFAbr','lLO','rLO'};
hemisphere = {'L'}; %, 'L', 'L', 'R'};

% Get the ROIs (actually just the spheres)
for iSub = 1:length(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];

    % for each region this subject has
    for iReg = 1:1 % length(mni{1}(:,1))

        % if the region is defined (VWFAbr for all, but in some is not present)
        % but we even need the VWFAbr ROI?
        if not(isnan(mni{iSub}(iReg, :)))

            % option 1: draw a sphere bounded by the mask
            % Sphere data: get the center, radius is fixed
            ROI_center = mni{iSub}(iReg, :);
            opt.radius = 10;

            % option 2: intersect localizer mask and neurosynth mask

            % Get the contrast to look for 
            switch iReg
                case 1, contrastToUse = 'frenchGtScrambled';
                otherwise, contrastToUse = 'drawingGtScrambled';
            end
            % Get the full name
            masksInDir = dir(fullfile(opt.dir.stats, subName, 'task-visualLocalizer_space-IXI549Space_FWHM-6', ...
                [subName, '_task-visualLocalizer_space-IXI549Space_desc-', contrastToUse ,'_*_p-0pt001_k-0_MC-none_mask.nii']));
            % Get the full path: folder + name
            localizerMask = fullfile(masksInDir.folder, masksInDir.name);

            % Get the reference, a random beta
            betaReference = fullfile(opt.dir.stats, subName, 'task-wordsDecoding_space-IXI549Space_FWHM-2', 'beta_0001.nii');

            % Get the neurosynth masks: loop through the options and make a lot of ROIs
            for iNS = 1:size(neurosynthMasks,2)
    
                % Get mask - already resliced to 
                nsMask = fullfile('masks', neurosynthMasks{iNS});

                % Open the image to make sure it's uint8
                % otherwise it's not read properly
                temp = load_nii(nsMask);
                temp.img = cast(temp.img, 'uint8');
                save_nii(temp, nsMask);

                % specify the sphere characteristics for each of them
                sphereParams = struct;
                sphereParams.location = ROI_center;
                sphereParams.radius = opt.radius;

                % specify the objects to pass: mask + sphere or mask + mask
                specSphere = struct('mask1', nsMask, ...
                                    'mask2', sphereParams);

                specMasks = struct('mask1', nsMask, ...
                                    'mask2', localizerMask);
    
                % specify the path for each subject
                outputPath = [opt.dir.rois,'/',subName];

                % make the mask - sphere 
                [sphereMask, sphereName] = createRoi('intersection', specSphere, betaReference, outputPath, opt.saveROI);

                % make the mask - masks
                [intersectedMask, intersectedName] = roi_createMasksOverlap(specMasks, betaReference, outputPath, opt.saveROI);
             
                % Rename .json and .nii files of both masks to have 
                % something more readable
                
                % Original names
                intersectedJustTheName = intersectedName(1:end-4);
                sphereJustTheName = sphereName(1:end-4);

                % New names
                intersectedNewName = fullfile(opt.dir.rois, subName, [subName, '_space-MNI_atlas-neurosynth_thresh-' ...
                                        num2str(iNS+4) '_method-masksIntersection_label-VWFAfr_mask']);
                sphereNewName = fullfile(opt.dir.rois, subName, [subName, '_space-MNI_atlas-neurosynth_thresh-' ...
                                    num2str(iNS+4) '_method-maskSphere_label-VWFAfr_mask']);

                % Rename sphere
                movefile(sphereName, [sphereNewName,'.nii'],'f')
                movefile([sphereJustTheName,'.json'], [sphereNewName,'.json'],'f')

                % Rename intersection
                movefile(intersectedName, [intersectedNewName,'.nii'],'f')
                movefile([intersectedJustTheName,'.json'], [intersectedNewName,'.json'],'f')

                % reslice the masks
                sphereMask = resliceRoiImages(betaReference, [sphereNewName, '.nii']);
                intersectedMask = resliceRoiImages(betaReference, [intersectedNewName, '.nii']);

            end

        end
    end
end

%% Method #7 - Anatomical masks: FG (BA37) + IT (BA20)
%
% BA masks from wfu_pickatlas 
% merged into one using custom script. 
% If deleted, run: 
% > roi_mergeMasks('masks/BA20.nii','masks/BA37.nii','masks/merged_BA20BA37')

% Get individual peaks
mni = roi_getMNIcoords(opt.subjects);

roiNames = opt.roiList;
save('ROI-trial_indivudal-coords.mat', 'roiNames', 'mni');

% which name to assign: 1 to 6 based on the region
regName = {'VWFAfr'}; %,'VWFAbr','lLO','rLO'};
hemisphere = {'L'}; %, 'L', 'L', 'R'};

% Get the ROIs (actually just the spheres)
for iSub = 1:length(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];

    % for each region this subject has
    for iReg = 1:1 % length(mni{1}(:,1))

        % if the region is defined (VWFAbr for all, but in some is not present)
        % but we even need the VWFAbr ROI?
        if not(isnan(mni{iSub}(iReg, :)))

            % option 1: draw a sphere bounded by the mask
            % Sphere data: get the center, radius is fixed
            ROI_center = mni{iSub}(iReg, :);
            opt.radius = 10;

            % option 2: intersect localizer mask and broadmann mask
            % Get the contrast to look for 
            switch iReg
                case 1, contrastToUse = 'frenchGtScrambled'; % VWFAfr/br
                otherwise, contrastToUse = 'drawingGtScrambled'; % r/lLO
            end
            % Get the full name
            masksInDir = dir(fullfile(opt.dir.stats, subName, 'task-visualLocalizer_space-IXI549Space_FWHM-6', ...
                [subName, '_task-visualLocalizer_space-IXI549Space_desc-', contrastToUse ,'_*_p-0pt001_k-0_MC-none_mask.nii']));
            % Get the full path: folder + name
            localizerMask = fullfile(masksInDir.folder, masksInDir.name);

            % Get the reference, a random beta
            betaReference = fullfile(opt.dir.stats, subName, 'task-wordsDecoding_space-IXI549Space_FWHM-2', 'beta_0001.nii');

            % Get BA mask
            baMask = 'masks/rBA37.nii';

            % Open the image to make sure it's uint8
            % otherwise it's not read properly
            temp = load_nii(baMask);
            temp.img = cast(temp.img, 'uint8');
            save_nii(temp, baMask);

            % specify the sphere characteristics for each of them
            sphereParams = struct;
            sphereParams.location = ROI_center;
            sphereParams.radius = opt.radius;

            % specify the objects to pass: mask + sphere or mask + mask
            specSphere = struct('mask1', baMask, ...
                                'mask2', sphereParams);

            specMasks = struct('mask1', baMask, ...
                               'mask2', localizerMask);

            % specify the path for each subject
            outputPath = [opt.dir.rois,'/',subName];

            % make the mask - sphere
            [sphereMask, sphereName] = createRoi('intersection', specSphere, betaReference, outputPath, opt.saveROI);

            % make the mask - masks
            [intersectedMask, intersectedName] = roi_createMasksOverlap(specMasks, betaReference, outputPath, opt.saveROI);

            % Rename .json and .nii files of both masks to have
            % something more readable

            % Original names
            intersectedJustTheName = intersectedName(1:end-4);
            sphereJustTheName = sphereName(1:end-4);

            % New names
            intersectedNewName = fullfile(opt.dir.rois, subName, ...
                [subName, '_space-MNI_atlas-Broadmann_method-masksIntersection_label-VWFAfr_mask']);
            sphereNewName = fullfile(opt.dir.rois, subName, ...
                [subName, '_space-MNI_atlas-Broadmann_method-maskSphere_label-VWFAfr_mask']);

            % Rename sphere
            movefile(sphereName, [sphereNewName,'.nii'],'f')
            movefile([sphereJustTheName,'.json'], [sphereNewName,'.json'],'f')

            % Rename intersection
            movefile(intersectedName, [intersectedNewName,'.nii'],'f')
            movefile([intersectedJustTheName,'.json'], [intersectedNewName,'.json'],'f')

            % reslice the masks
            sphereMask = resliceRoiImages(betaReference, [sphereNewName, '.nii']);
            intersectedMask = resliceRoiImages(betaReference, [intersectedNewName, '.nii']);

        end
    end
end


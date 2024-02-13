%% Create ROIs based on peak coordinates from subjects' localizers
%
% From Fedorenko's parcels (Fedorenko et al. 2010).
% Extract the language ROIs, reslice them, and apply them to individual data
%
% Outputs:
% - ROIs (duh) in the 'code/rois/masks'
% - ROIs in the folder specified in roi_option

% If not done previously, extract language ROIs from Fedorenko's parcels
if isempty(dir('masks/fedorenko_parcels/r*'))

    getRoiFromParcels(opt);
end

%% Overlap the fROI to each subject VWFA contrast 
% From each subject's constrasts [FW > SFW] and [BW > SBW], overlap these 
% contrasts with each of Fedorenko's fROIs.
% Outputs a report 

% Get the resliced masks in the folder (created by code above)
fedorenkoMasks = dir('masks/fedorenko_parcels/r*');

% Set an accepted nb of voxels 
nbVoxels = opt.numLanguageVoxels;

% Initialize report 
report = {'subject','script','hemi','area','done','voxels','enough'};


for iSub = 1:length(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];

    % Get a reference image for this sub
    dataImage = fullfile(opt.dir.stats, subName, ...
                         ['task-visualLocalizer_space-',opt.space{1},'_FWHM-6_node-localizerGLM'], 'beta_0001.nii');

    % Get the contrasts: [FW > SFW] 
    % Only works on FR contrast. Assumption is to find BR effects in
    % FR-defined areas
    subConFR = dir(fullfile(opt.dir.stats, subName, ['task-visualLocalizer_space-',opt.space{1},'_FWHM-6_node-localizerGLM'], ...
                   ['sub-*_space-',opt.space{1},'_desc-f*pt05*_mask.nii']));

    % Join them to make looping easier
    subContrasts = vertcat(subConFR);

    thisContrast = fullfile(subContrasts.folder, subContrasts.name);
    
    % Get the type of stimuli of the contrast, to be specified later in
    % the new ROI name
    conName = strsplit(thisContrast,{'desc-','Gt'});
    con = conName{2};

    % From this contrast, cut out all the necessary fROIs
    for iFed = 1:size(fedorenkoMasks,1)
        
        thisROI = fullfile(fedorenkoMasks(iFed).folder, fedorenkoMasks(iFed).name);

        % Reslice on each sub, necessary for PPI.
        % Only one time, it's not necessary to do it for each contrast
        unslicedRoi = fullfile(fedorenkoMasks(iFed).folder, fedorenkoMasks(iFed).name(2:end));
        resliceOnParticipant(unslicedRoi, opt, subName);

        % Get name and hemisphere
        % (highlighly sensitive to path and folders)
        strName = strsplit(thisROI, {'_','-'});
        hemi = strName{5};
        reg = strName{11};

        % Load the mask and cast it as uint8, 
        % otherwise it's not recognized as a binary mask
        recastROI = load_nii(thisROI);
        recastROI.img = cast(recastROI.img, 'uint8');
        save_nii(recastROI, thisROI);
        
        % specify the objects to pass: mask + mask
        specMasks = struct('mask1', thisROI, ...
                           'mask2', thisContrast);

        % specify the path for each subject
        outputPath = fullfile(opt.dir.rois, subName);

        % Reset fROI name. In case it stays empty, it means the ROI was
        % not created, thus there is no name change to do
        froiName = [];

        % Intersection of localizer spmT and mask (TBD)
        [froiMask, froiName] = roi_createMasksOverlap(specMasks, dataImage, outputPath, opt.saveROI);

        % Only rename ROI if an ROI was indeed created in the previous
        % line, otehrwise try with the next one
        if ~isempty(froiName)
            % Rename .json and .nii files of both masks to have more readable names
            % Remove file extension from name
            froiJustName = froiName(1:end-4);
            % New names
            froiNewName = fullfile(opt.dir.rois, subName, [subName,'_hemi-',hemi, ...
                                   '_space-',opt.space{1},'_atlas-language_contrast-',con,'_label-',reg,'_mask']);
            % Rename intersection
            movefile(froiName, [froiNewName,'.nii'],'f')
            movefile([froiJustName,'.json'], [froiNewName,'.json'],'f')
            % reslice the masks
            intersectedMask = resliceRoiImages(dataImage, [froiNewName, '.nii']);
        end

        % Add information to the report
        % - whether ROI was created or skipped (empty)
        if isempty(froiName) 
            done = 'skipped';
        else
            done = 'created';
        end

        % - were there enough (more than the threhold number of) voxels? 
        if froiMask.roi.size >= nbVoxels
            enough = true;
        else
            enough = false;
        end

        % - how many voxels are there in the new ROI?
        voxels = froiMask.roi.size;

        % Add subject, contrast, hemisphere, region, and details to report
        report = vertcat(report, ...
                         {subName, con, hemi, reg, done, voxels, enough});

    end
end

% Save report
writecell(report,['reports/roi_reports_languageROIs_voxThres-' num2str(nbVoxels) '_' date '.txt']);


%% FUNCTION TO EXTRACT ALL FEDORENKO'S PARCELS

function getRoiFromParcels(opt)
% Extract ROIs from atlas is taken and adapted from the
% extractRoiFromAtlas function of CPP_ROI

% load look-up table and get ROI names
fedLut = spm_load('masks/fedorenko_parcels/LUT.csv');
roiNames = fedLut.ROI';

% load parcels .nii file and get the volume data
fedHdr = spm_vol('masks/fedorenko_parcels/allParcels_language_SN220.nii');
fedVol = spm_read_vols(fedHdr);

% There are 12 ROIs (see lut.csv for a list of them)
% Loop through them and extracct them all
for r = 1:numel(roiNames)

    label = fedLut.label(r);
    labelStruct = struct('ROI', roiNames{r}, ...
                         'label', label);
    
    % create 'negative' vol that is all zeros, set as one only the voxels
    % corresponding with the current ROI label
    outputVol = false(size(fedVol));
    outputVol(fedVol == labelStruct.label) = true;

    % Personal addition, adhere to the common naming for ROIs
    % Split the name of the parcel in two: 
    % - hemisphere goes after 'hemi' attribute,
    % - area goeas after 'label'
    strName = strsplit(roiNames{r}, '_');
    hemi = strName{1}(1);
    reg = strName{2};
    roiFilename = ['hemi-',hemi,'_space-',opt.space{1},'_atlas-language_label-',reg,'_mask.nii'];

    % Assign the new filename to the ROI and save it
    fedHdr.fname = spm_file(fedHdr.fname, 'filename', roiFilename);
    spm_write_vol(fedHdr, outputVol);

    % reslice them to be in the same space as our beta maps
    % to do so, we need a beta map of reference
    % It does not matter which sub now, when we apply them to single
    % subjecct they'll be resliced with their distorsions
    dataImage = fullfile(opt.dir.stats, 'sub-006', ...
                         ['task-visualLocalizer_space-',opt.space{1},'_FWHM-6_node-localizerGLM'], 'beta_0001.nii');

    resliceRoiImages(dataImage, fullfile('masks/fedorenko_parcels', roiFilename));

end
end

%% FUNCTION TO RESLICE ROI BASED ON PARTICIPANT'S SPACE 
function resliceOnParticipant(roi, opt, subName)

    % copy image in sub roi folder
    roiParseName = strsplit(roi, 'fedorenko_parcels/');
    copyfile(roi, fullfile(opt.dir.rois, subName, [subName '_' roiParseName{2}]), 'f');
    
    % get new name
    copiedRoi = fullfile(opt.dir.rois, subName, [subName '_' roiParseName{2}]);

    % Open the mask, cast as binary, close it 
    recastROI = load_nii(copiedRoi);
    recastROI.img = cast(recastROI.img, 'uint8');
    save_nii(recastROI, copiedRoi);

    % get reference
    dataImage = fullfile(opt.dir.stats, subName, ...
                         ['task-visualLocalizer_space-',opt.space{1},'_FWHM-6_node-localizerGLM'], 'beta_0001.nii');

    % relisce roi based on the specific sub
    resliceRoiImages(dataImage, copiedRoi);

    % case resliced ROI as binary 
    reslicedRoi = fullfile(opt.dir.rois, subName, ['r' subName '_' roiParseName{2}]);
    % Open the mask, cast as binary, close it 
    recastROI = load_nii(reslicedRoi);
    recastROI.img = cast(recastROI.img, 'uint8');
    save_nii(recastROI, reslicedRoi);

end



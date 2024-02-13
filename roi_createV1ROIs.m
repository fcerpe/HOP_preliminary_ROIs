%% Create ROIs for V1 area 
%
% Extract ROIs for eacch participant and subject, reslice them, and apply 
% them to individual data
% From anatomy toolbox and extracted V1 mask.

% Define the atlas
% visfatlas was previously implemented
atlas = 'jubrain';
visfatlasV1Masks = dir('masks/anatomy toolbox/*V1*');

% Extract mask, actually the only one  
v1mask = fullfile(visfatlasV1Masks(1).folder, visfatlasV1Masks(1).name);

% Iterate through subjects to reslice the original mask on each subject's
% space and save it as a personal ROI
for iSub = 1:numel(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];

    % Get a reference image for this sub
    dataImage = fullfile(opt.dir.stats, subName, ...
                         ['task-visualLocalizer_space-',opt.space{1},'_FWHM-6_node-localizerGLM'], 'beta_0001.nii');

    % Get name and use it to compose new path
    newMaskPath = fullfile(opt.dir.rois, subName, [subName '_' visfatlasV1Masks(1).name]);
    copyfile(v1mask, newMaskPath, 'f');

    % Open the mask and cast it as binary
    recastROI = load_nii(newMaskPath);
    recastROI.img = cast(recastROI.img, 'uint8');
    save_nii(recastROI, newMaskPath);

    % Relisce roi based on the specific sub
    resliceRoiImages(dataImage, newMaskPath);

    % Get resliced ROI as binary 
    reslicedRoi = fullfile(opt.dir.rois, subName, ['r' subName '_' visfatlasV1Masks(1).name]);

    % Open the mask and cast it as binary
    recastROI = load_nii(reslicedRoi);
    recastROI.img = cast(recastROI.img, 'uint8');
    save_nii(recastROI, reslicedRoi);

    % Get the contrasts: [FW + SFW > nothing] 
    % And the full path
    subCon = dir(fullfile(opt.dir.stats, subName, ['task-visualLocalizer_space-',opt.space{1},'_FWHM-6_node-localizerGLM'], ...
                   'sub-*_desc-allF*pt05*_mask.nii'));
    thisContrast = fullfile(subCon.folder, subCon.name);
    
    % Get the type of stimuli of the contrast, to be specified later in
    % the new ROI name
    conName = strsplit(thisContrast,{'desc-','Gt'});
    con = conName{2};

    % Get name and hemisphere
    strName = strsplit(v1mask, {'_','-'});
    hemi = strName{4};
    reg = strName{10};

    % Load the mask and cast it as uint8
    recastROI = load_nii(v1mask);
    recastROI.img = cast(recastROI.img, 'uint8');
    save_nii(recastROI, v1mask);

    % specify the objects to pass: mask + mask
    specMasks = struct('mask1', reslicedRoi, ...
                       'mask2', thisContrast);

    % specify the path for each subject
    outputPath = fullfile(opt.dir.rois, subName);

    % Reset fROI name. In case it stays empty, it means the ROI was
    % not created, thus there is no name change to do
    froiName = [];

    % Intersection of localizer spmT and mask (TBD)
    [froiMask, froiName] = roi_createMasksOverlap(specMasks, dataImage, outputPath, opt.saveROI);

    % Only rename ROI if an ROI was indeed created in the previous line, 
    % otherwise try with the next one
    if ~isempty(froiName)

        % Rename both .json and .nii files of both masks to have more readable names
        froiJustName = froiName(1:end-4);

        % Custom names
        froiNewName = fullfile(opt.dir.rois, subName, [subName,'_hemi-',hemi, ...
                               '_space-',opt.space{1},'_atlas-',atlas,'_contrast-',con,'_label-',reg,'_mask']);

        % Rename intersection
        movefile(froiName, [froiNewName,'.nii'],'f')
        movefile([froiJustName,'.json'], [froiNewName,'.json'],'f')

        % Reslice the mask
        intersectedMask = resliceRoiImages(dataImage, [froiNewName, '.nii']);
    end

end




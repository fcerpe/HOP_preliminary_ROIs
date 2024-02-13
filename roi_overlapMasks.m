%% Calculate VWFA overlap between-subjects
%
% Take the ROIs for each subject and sum them
% Display the overlap as % of shared space

% Create a struct to store all the rois. Left empty for now and initalized
% while loading the first roi
vwfas = [];

for iSub = 1:length(opt.subjects)

    % Get subject number
    subName = ['sub-', num2str(opt.subjects{iSub})];

    fprintf(['Adding vwfa for ' subName '\n']);

    % Get the subject's VWFA mask
    brainMask = fullfile(opt.dir.rois, subName, ...
        ['r' subName '_hemi-L_space-MNI_atlas-neurosynth_method-expansionIntersection_label-VWFAfr_mask.nii']);

    % Load the mask
    % USe load_nii instead of our custom function because we deal with
    % multiple subjects at the same time. A more simple measure appears to
    % be more solid
    vwfaMask = load_nii(brainMask);

    % In case it's the first mask, use it as base to summ all the others
    % Otherwise, add this .img to the previous masks
    if isempty(vwfas)
        vwfas = struct;
        allImg = vwfaMask.img;
        totalMasks = 1;
    else
        allImg = allImg + vwfaMask.img;
        totalMasks = totalMasks + 1;
    end

    % Add it to the others, for storage
    eval(['vwfas.sub' num2str(opt.subjects{iSub}) ' = vwfaMask;']);

end

% Create a custom mask with the overlap
vwfaMaskName = fullfile(opt.dir.rois, 'overlap_hemi-L_space-MNI_label-VWFAfr.nii');
vwfaMask.img = allImg;

save_nii(vwfaMask, vwfaMaskName);

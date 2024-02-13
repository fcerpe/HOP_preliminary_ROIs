% Process report

d = '18-Sep-2023_voxThres-50';

% Don't procastiante, get today's report!
repFiles = dir(['languageRoiReport_' d '*.txt']);

% load last file, the most recent
report = readtable(fullfile(repFiles(end).folder, repFiles(end).name));

% Get only the relevant ROIs
% - more than 50 voxels
% - skim some columns: no script or created necessary
enoughVoxels = report(report.enough == 1, [1 3 4 6]);

% Get the number of participants processed, to obtain a measure of the
% consensus in the activations
nbParticipants = numel(unique(report.subject));

% Get a list of the different areas
% load look-up table and get ROI names
hemis = {'L','R'};
areas = {'IFGorb','IFG','MFG','AntTemp','PosTemp','AngG'};
recap = table('Size',[12 4],'VariableTypes',{'string','categorical','double','double'},'VariableNames',{'hemi','area','nbROIs','percentage'});

% Scroll through all hemispheres and areas to lookup how many participants
% have sufficiently large clusters in the ROI.

rec = 1;

for h = 1:numel(hemis)
    thisHemi = hemis{h};
    eval(['rois_' thisHemi ' = enoughVoxels(find(ismember(enoughVoxels.hemi, ''' thisHemi ''')), :);']);

    for a = 1:numel(areas)
        thisArea = areas{a};
        eval(['rois_' thisHemi '_' thisArea ' = rois_' thisHemi '(find(ismember(rois_' thisHemi '.area, ''' thisArea ''')), :);']);
        eval(['nbROIs = size(rois_' thisHemi '_' thisArea ',1);']);

        % get who has those activations: experts or controls
        eval(['participants = rois_' thisHemi '_' thisArea '.subject;']);

        nbExp = 0;
        nbCtr = 0;

        for p = 1:numel(participants)
            switch participants{p}(end-1:end)
                case {'06', '07', '08', '09', '12', '13'}
                    nbExp = nbExp + 1;
                otherwise
                    nbCtr = nbCtr + 1;
            end
        end

        % add data to recap table
        recap.hemi(rec) = thisHemi;
        recap.area(rec) = [thisHemi, ' ', thisArea];
        recap.nbROIs(rec) = nbROIs;
        recap.nbExp(rec) = nbExp;
        recap.nbCtr(rec) = nbCtr;
        recap.percentage(rec) = nbROIs / nbParticipants;
        rec = rec +1;

    end
end

%% Plot results and save report image

bar(recap.area, [recap.nbExp, recap.nbCtr], 'stacked')

% individual numbers
yyaxis left
ylim([0 nbParticipants])
yticks(0:nbParticipants)
ylabel('nb of subjects')

% total percentages
yyaxis right
ylim([0 1])
ax = gca;
ax.YColor = 'b';
ylabel('total percentage of subjects')
yticks(0:0.1:1)
yline(0.8,'--r')

% legend
legend('experts','controls')

% save fig
saveas(gcf,['languageRoiReport_' d '.png'])


% Print on screen decision about which areas to use.
% IMPORTANT: this is purely visual. 
%            In ppi_option and mvpa_option, the areas on which to perform
%            analyses and subject pool are coded manually
% TO-DO (14/09/2023)
% - extend communication between scripts and make this decision known to
%   other analyses (maybe by loading a txt file)

fprintf(['\n\nThe following areas are present in more than 80%% of the subjects:\n\n']);
for ar = 1:size(recap,1)
    if recap.percentage(ar) >= 0.8
        fprintf(['\t- ' char(recap.area(ar)) '\n']);
    end
end
fprintf(['\nWill use the ROIs specified in the list to probe for higher level language effects\n' ...
         '(using PPI and MVPA)\n']);


%% Create overlap of the selected areas 
% from roi_splitVWFA: calculate VWFA overlap between-subjects

% Create a struct to store all the rois. Left empty for now and initalized
% while loading the first roi

fprintf(['\n\n Creating overlap masks\n']);


selectedAreas = recap(recap.percentage > .7, :);

leftPosTemp = [];
rightPosTemp = [];

for sa = 1:size(selectedAreas,1)
    areaToGet = char(selectedAreas.area(sa));
    eval(['currArea = rois_' selectedAreas.hemi{sa} '_' areaToGet(3:end) ';']);


    eval([selectedAreas.hemi{sa} '_' areaToGet(3:end) 's = [];']);
    tempRegistry = [];

    for iSub = 1:length(currArea.subject)

        % Get subject number
        subName = currArea.subject{iSub};

        fprintf(['Adding ' areaToGet ' for ' subName '\n']);

        % Get the subject's VWFA mask
        brainMask = fullfile(opt.dir.rois, subName, ...
            ['r' subName '_hemi-' selectedAreas.hemi{sa} '_space-MNI_atlas-fedorenko_contrast-french_label-PosTemp_mask.nii']);

        % Load the mask
        % USe load_nii instead of our custom function because we deal with
        % multiple subjects at the same time. A more simple measure appears to
        % be more solid
        currMask = load_nii(brainMask);

        % In case it's the first mask, use it as base to summ all the others
        % Otherwise, add this .img to the previous masks
        if isempty(tempRegistry)
            tempRegistry = struct;
            allImg = currMask.img;
            totalMasks = 1;
        else
            allImg = allImg + currMask.img;
            totalMasks = totalMasks + 1;
        end

        % Add it to the others, for storage
        eval(['tempRegistry.sub' num2str(opt.subjects{iSub}) ' = currMask;']);

    end

    % Create a custom mask with the overlap
    currMaskName = fullfile(opt.dir.rois, ['overlap_hemi-' selectedAreas.hemi{sa} '_space-MNI_label-' areaToGet(3:end) '.nii']);
    currMask.img = allImg;

    save_nii(currMask, currMaskName);

    eval([selectedAreas.hemi{sa} '_' areaToGet(3:end) 's = tempRegistry;']);

end







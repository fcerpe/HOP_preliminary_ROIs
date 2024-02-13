%% Solve overlapping ROIS
% In the case that two ROIs have voxels in common, assign them to the
% closest peak

% Prepare report 
opt.peaks = struct('sub','001','method','center','area','none','coords',[0,0,0]);

iPeak = 1;

% This script only checks VWFA defined by [FW > SFW] contrast and lLO, the
% only two areas with possible overlap
opt.overlapsToCheck = {'VWFAfr-lLO'};

for iSub = 1:numel(opt.subjects)

    % Skip subjects 23 and 24, no lLO activation
    if ~strcmp(opt.subjects{iSub}, '023') && ~strcmp(opt.subjects{iSub}, '024')

        subName = ['sub-', opt.subjects{iSub}];
    
        % Notify the user 
        fprintf(['Analysing ' subName '\n']);
    
        % Load ROIs: 1-VWFA, 2-lLO

        % Get names to overwrite them later
        roi1name = ['r',subName,'_hemi-L_space-',opt.space{1},'_atlas-neurosynth_method-expansionIntersection_label-VWFAfr_mask.nii'];
        roi2name = ['r',subName,'_hemi-L_space-',opt.space{1},'_atlas-neurosynth_method-expansionIntersection_label-lLO_mask.nii']; 

        roi1 = load_nii(fullfile(opt.dir.rois, subName, roi1name));
        roi2 = load_nii(fullfile(opt.dir.rois, subName, roi2name));
    
        % Load contrasts, to resolve cases of equal distance from peaks.
        % In the case of equal distances AND different contrasts (VWFA - lLO),
        % assign the voxel to the contrast with the most activity
        spmT_words = dir(fullfile(opt.dir.stats, subName, ['task-visualLocalizer_space-',opt.space{1},'_FWHM-6_node-localizerGLM'], ...
                         [subName, '_*_desc-frenchGtScrambled_p-0pt001_*_mask.nii']));
        spmT_objects = dir(fullfile(opt.dir.stats, subName, ['task-visualLocalizer_space-',opt.space{1},'_FWHM-6_node-localizerGLM'], ...
                           [subName, '_*_desc-drawingGtScrambled_p-0pt001_*_mask.nii']));
    
        conWords = fullfile(spmT_words.folder, spmT_words.name);
        conObjects = fullfile(spmT_objects.folder, spmT_objects.name);
    
        %% Preliminary calculations
    
        % Get the subject's peaks
        mni = roi_getMNIcoords(opt.subjects);
    
        % Sum the masks: 
        % value of 2 indicate activation in both areas
        overlapMask = roi1.img + roi2.img;
        
        % VWFA peak
        peakWords = mni{iSub}(1,:);
    
        % lLO peak
        peakObjects = mni{iSub}(3,:);
    
        % get areas and peaks positions for both ROIs
        roi1area = find(roi1.img > 0);
        roi2area = find(roi2.img > 0);
    
        % take the indexes of the overlapping cells (voxels)
        overlap = find(overlapMask == 2);
        [overlap(:,1), overlap(:,2), overlap(:,3)] = ind2sub([61,73,61],overlap);
    
        if not(isempty(overlap))
    
            % Notify the user
            fprintf(['Overlap between VWFA and lLO: ' num2str(size(overlap,1)) ' voxels\n']);
    
            % Assign voxels to either ROI
            % go through all the overlap voxels, look at which peak is closest 
            % and assign the voxel to that ROI (cancel it from the other ROI)
            % If they are at the same distance, look at univariate activation 
            % and assign the voxel to the highest bidder
    
            % nbIstances indicates how many ties there are, for reference
            nbIstances = 0;
    
            for ov = 1:size(overlap,1)
                currVoxelCoords = overlap(ov,:);
    
                % distance from the peaks
                distanceFromWords = abs(currVoxelCoords - peakWords);
                distanceFromWords = sum(distanceFromWords);
    
                distanceFromObjects = abs(currVoxelCoords - peakObjects);
                distanceFromObjects = sum(distanceFromObjects);
    
                % Assigning to one area means removing the corresponding value
                % from the opponent
                if distanceFromWords < distanceFromObjects
                    roi2.img(currVoxelCoords(1),currVoxelCoords(2),currVoxelCoords(3)) = 0;
    
                elseif distanceFromWords > distanceFromObjects
                    roi1.img(currVoxelCoords(1),currVoxelCoords(2),currVoxelCoords(3)) = 0;
    
                else
                    nbIstances = nbIstances +1;
    
                    % We can compare activations between two different contrasts
                    % decide based on univariate response 
    
                    % If in that voxel, activation for words is higher than
                    % activation for objects
                    if conWords.img(currVoxelCoords(1), currVoxelCoords(2), currVoxelCoords(3)) >  ...
                            conObjects.img(currVoxelCoords(1), currVoxelCoords(2), currVoxelCoords(3))
    
                        roi2.img(currVoxelCoords(1),currVoxelCoords(2),currVoxelCoords(3)) = 0;
    
                    else
                        roi1.img(currVoxelCoords(1),currVoxelCoords(2),currVoxelCoords(3)) = 0;
                    end
                end
            end
    
            % Notify the user
            fprintf(['Re-assigned ' num2str(size(overlap,1)) ' voxels, of which ' num2str(nbIstances) ' based on activation\n']);
    
            % Post-transformation check
            % sum the masks again
            overlapMask = roi1.img + roi2.img;
            roi1area = find(roi1.img > 0);
            roi2area = find(roi2.img > 0);
    
            % take the indexes of the overlapping cells (voxels)
            overlap = find(overlapMask == 2);
            [overlap(:,1), overlap(:,2), overlap(:,3)] = ind2sub([61,73,61],overlap);
    
            % Notify the user
            fprintf(['There are now ' num2str(size(overlap,1)) ' overlapping voxels between VWFA and lLO\n\n\n']);
    
            % Save new ROIs
            save_nii(roi1,fullfile(opt.dir.rois, subName, roi1name));
            save_nii(roi2,fullfile(opt.dir.rois, subName, roi2name));
    
        else
            % Notify the user
            fprintf('No overlap between VWFA and lLO\n\n\n');
    
        end
    end

end





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
        
        % Load T constrast 
        [tConImg, tStruct] = loadTContrastImage(model, currROI.contrast, subjectFolder);

        
        % Subfunction to create expansion - not yet implemented
        % Will just show workflow for now, names may not match
        
        % call subfunction to expand
        %   get details about tmap (thersholds, map itself)
        %   while not enough voxels and enough thresholds
        %       apply threshold to tmap
        %       call subfunction to expand around a tmap
        %           get details about sphere (center, radius, stesp, maximum expansion)
        %           while not enough voxels and enough space to expand
        %               call subfunction to make sphere
        %               call subfunction to intersect sphere and tmap
        %               compare number of voxels
        %               if ok, retrun success = true and roi
        %           if not enough space, return success = false
        %       if success, return roi
        %       if not success, go to new threshold
        %   if never success, throw error
        
        % Save the resulting ROI to files(s)
        % saveROI(roi, finalPThresh, currROI.task, currROI.contrast, outputFolder);

    end
end

end


%% Subfunctions

% make expansion
function roi = expandOnTContrast(opt, m, currROI, tConImg)

% List all the possible threshold steps that we can use 
% adjustmentPThreshs = [0.001, 0.01, 0.05];

% control for nbVoxels
% currNbVoxels = 0;

% Loops until certain conditions are met:
% - the target number of voxels is reached
% - the threshold is too relaxed
% while (currNbVoxels < m.targetNbVoxels) && (iThr <= length(adjustmentPThreshs))

    % threshold image
    % [thresholdedImage, numCells] = applyThreshold(tConImg, finalasPThresh, erdf);

    % Call subfunction to work on one threshold
    % [success, overlap] = expandOnThreshold(opt, m, currROI, thresholdedImage)



% end


%     % Notify the user
%     fprintf(1, '\n radius: %0.2f mm; roi size: %i voxels', radius, currNbVoxels);
% 



end

function [roi, success] = expandOnThreshold(opt, m , currROI, tConImg)

% Get specifics of radius 
% - maximum
% - steps based on voxel size
% - index

% while (currNbVoxels < m.targetNbVoxels) && (mask.roi.radius <= maxRadius)

    % Call subfunction to make a sphere
    % mask = maroi_sphere(struct('centre', cetner, 'radius', radius, ...
    %                            'label', 'expansion', 'binarize', 1, 'roithresh', 0.5));

    % Call subfunction to intersect sphere and mask

    % compare number of voxels with target

    % Check that everything is in order
    % - size is indeed expanding. If not, it means that the sphere has 
    %   engulfed the whole cluster of voxels 
    % - there are still voxels to include in the ROI. If not, it means that
    %   the threshold is too strict
    % - the radius is larger than the mask itself. If so, there are no more
    %   voxels in the brain
    % if (mask.roi.size == previousSize) || (mask.roi.size > sphere.maxNbVoxels) || (mask.roi.radius > maxRadius)   
    %   iThr = iThr + 1;
    % end

end
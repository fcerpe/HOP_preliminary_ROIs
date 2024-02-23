%% HOPLAB fMRI pipeline - preliminary ROI creation 
%
% Main script for the creation of ROIs. From the options inputed in 'roi_option()', 
% will follow different methods and produce preliminary rois. 
%
% [Andrea's way of writing doc] 


%% Load options and necessary toolboxes 

% nice and clean
clear
clc

% Initialize marsbar, spm, etc
% (all taken care of by bidspm)
addpath 'lib/bidspm'
bidspm;

% Load options from function
[opt, roiMethods] = hop_option();


%% extract preliminray ROIs

% Follows parameters defined by the user in 'hop_option()'.
% Several methods and their specifics are stored in the struct 'roi_methods'

for iM = 1:numel(roiMethods)

    % Notify the user
    fprintf(['\n\n Processing method #', num2str(iM), '\n\n']);

    % Extract ROIs ased on the method indicated
    switch roiMethods(iM).method

        case 'intersection'
            % Intersect a spmT map with either a sphere or another mask
            hop_createROI_intersection(opt, roiMethods(iM));

        case 'sphere'
            % Create a sphere around specified coordinates
            ac_makeSphereROI(opt, roiMethods(iM)); % -> hop_createROI_sphere

        case 'atlas'
            % Take region from a specified area of an atlas
            hop_createROI_atlas(opt, roiMethods(iM));

        case 'expansion'
            % Expand a progressively large sphere from specified
            % coordinates and intersect it with a specified spmT map
            hop_createROI_expansion(opt, roiMethods(iM));

    end

end
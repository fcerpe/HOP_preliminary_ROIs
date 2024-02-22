%% HOPLAB fMRI pipeline - preliminary ROI creation 
%
% Main script for the creation of ROIs. From the options inputed in 'roi_option()', 
% will follow different methods and produce preliminary rois. 
%
% [Andrea's way of writing doc] 


%% Load options and necessary toolboxes 

% nice and clean
clear;
clc;

% Initialize marsbar, spm, etc
% (all taken care of by bidspm)
addpath '../lib/bidspm'
bidspm;

% Load options from function
[opt, roiMethods] = hop_option();


%% extract preliminray ROIs

% Follows parameters defined by the user in 'hop_option()'.
% Several methods and their specifics are stored in the struct 'roi_methods'

for iM = 1:numel(roi_methods)

    % Extract ROIs ased on the method indicated
    switch roi_methods(iM).method

        case 'intersection'
            % Intersect a spmT with either a sphere or another mask
            hop_createROI_intersection(opt);

        case 'sphere'
            % Create a sphere around specified coordinates
            hop_createROI_sphere(opt);

        case 'atlas'
            % Take region from a specified area of an atlas
            hop_createROI_atlas(opt);

        case 'anatomical'
            hop_createROI_anatomical(opt);

        case 'expansion'
            hop_createROI_expansion(opt);

    end

end
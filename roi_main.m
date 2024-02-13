%% VISual BRAille ROI creation 
%
% main script of the pipeline, extracts regions of interest (ROIs) for all
% the regions we are, well, interested in.
% Uses different scripts / methods based on the avialbility of localizer
% data. 
% The use of different methods is due to a limited localizer experiment
% IMPORTANT: if data is stored using datalad, you need to unlock the
% 'stats' and 'rois' folders first
% 
% ROIs and methods
% - VWFA, lLO, rLO: localzier contrast [FW > SFW], expansion intersection
%
% - V1: intersection between [FW + SFW > rest] and anatomy toolbox V1 area
%
% - l-PosTemp (and other language areas): intersection between localizer 
%   contrast [FW > SFW] and Fedorenko's atlas 


clear;
clc;

% Get path and init bidspm
addpath '../lib/bidspm'
bidspm;

% Get options - common to all the scripts
opt = roi_option();


%% Select the ROIs based on the functional localizers 
% Extract the following ROIs based contrasts coming from the functional localizers
% - VWFAfr [french words]  - [scrambled french words]
% - lLO    [line drawings] - [scrambled line drawings]
% - rLO    [line drawings] - [scrambled line drawings]
 
% Extract ROIs
roi_createROIs; 

% Resolve possible overlaps between neighbouring areas
roi_resolveOverlaps; 

% Attempt to split VWFA into posterior / anterior sub-areas. 
% Result is non-conclusive, but script kept for posterity
roi_splitVWFA;


%% Extract language ROIs 
% From Fedorenko et al.'s parcels
%
% Create subject level ROIs:
% - just the mask resliced to the participant's reference
% - intersections of mask and [FW > SFW]
%
% Additionally, generate a report about how many / which subjects present 
% activation in which parcels. 
% This report drives the decision of which areas to focus on in PPI
% (more in code/ppi)

roi_createLanguageROIs;


%% Extract V1 ROIs 
% Similar process to language ROIs, takes a V1 ROI form Anatomy Toolbox and 
% overlaps it to the the mask for [FW + SFW > rest] contrast

roi_createV1ROIs;


%% Create overlap masks
% For each ROI, create an overlap of each subject's mask

% roi_overlapMask; 




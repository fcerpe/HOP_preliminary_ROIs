function [opt, roiMethods] = hop_option()

%% HOPLAB ROI options 
%
% Stores the options created from the user to make new ROIs
%
% [Andrea's way of writing doc] 
%
% returns two structures:
%       
%       opt: options for the different ROIs and some method-general paramteres (TBD)
%       roiMethods: table-like structure that specifies all the methods and
%                   the main parameters required to run the selected
%                   extractions

% Avoid possible overlaps for multiple runs 
opt = [];
roiMethods = [];


%% Main parameters
% 
% these parameters are valid for all the methods and are stored in 'opt'

% Paths 
opt.dir.root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
opt.dir.raw = fullfile(opt.dir.root, 'inputs', 'raw');
opt.dir.derivatives = fullfile(opt.dir.root, 'outputs', 'derivatives');
opt.dir.preproc = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-preproc');
opt.dir.input = opt.dir.preproc;
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'cpp_spm-rois');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-stats');
opt.dir.cosmo = fullfile(opt.dir.root, 'outputs', 'derivatives', 'CoSMoMVPA');
opt.dir.jobs = fullfile(opt.dir.stats, 'jobs', opt.taskName);


% Space of the data 
% - IXI549Space: bidspm preprocessing
% - MNI152NLin2009cAsym: fmriprep preprocessing
% - individual: just in case it's needed
opt.space = 'MNI152NLin2009cAsym'; 


% Subjects on which to perform ROI extraction
% (default value, applied if not specified otherwise)
opt.subjects = {'006','007'}; 


% Task to analyze
opt.taskName = 'task';


%% Template methods 
% 
% Each method serves as a template
% Feel free to modify / copy / delete the following methods to your needs,
% adding all the subjects, areas, coords, contrasts needed


% METHOD #1 - Sphere of [radius]mm in [area]
roiMethods(1).method = 'sphere';

% on which subjects should we extract ROIs?
roiMethods(1).subjects = {'006','007'}; 

% what radii should the sphere(s) have?
roiMethods(1).radii = [10, 8, 6];

% which areas will we extract? 
roiMethods(1).roiNames = {'area1';
                          'area2'};

% from which coordinates? 
% (if lateralized, omit controlateral coords using NaNs)
roiMethods(1).coordsLeft =  [0 0 0; 
                             15 15 15];
roiMethods(1).coordsRight =  [0 0 0; 
                              15 15 15];



% METHOD #2 - Intersection between [contrast] and [roi]
roiMethods(2).method = 'intersection';

% where to get the ROI and the tmap?
roiMethods(2).roiPath = '../whatever';
roiMethods(2).tmapPath = '../whatever';

% which min. number of voxel can we accept? 
roiMethods(2).nbVoxel = 25;



% METHOD #3 - Atlas: extraction of [area] from [atlas]
roiMethods(3).method = 'atlas';

% Which atlas are we using?
% - select from list
% - give for granted that the masks are in order
roiMethods(3).atlas = 'Broadmann';
roiMethods(3).mask = 'BA8';
roiMethods(3).maskPath = '../whatever';




% METHOD #4 - Expansion around [area] in [contrast]
roiMethods(4).method = 'expansion';
roiMethods(4).nbVoxel = 100;



% description to add to folder name, to distinguish from GLM 
% Possibly obsolete
opt.desc = 'ROI';




%% Check that everything is in order

% Sphere - all parameters should be present

% Intersection - all parameters should be present

% Expansion - all parameters should be present

% Atlas - all parameters should be present and the atlas should be present


end


%% Support functions



outList = struct;

% Assign method
outList.method = inList.method;

% Based on the method, combine other parameters
switch inList.method

    case 'sphere'
        % make combinations

        % assign coordinates to corresponding area


end

end
















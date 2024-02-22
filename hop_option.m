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



%% Template methods 
%
% Feel free to modify / copy / delete the following methods to your needs


% Method #1 - Sphere


% Method #2 - Intersection


% Method #3 - Atlas 


% Method #4 - 



% suject to run in each group
opt.subjects = {'006','007','008','009','010','011','012','013','018','019','020','021','022','023','024','026','027','028'}; 
% Avaiable participants: 
% '006','007','008','009','010','011','012','013','018','019','020','021','022','023','024','026','027','028'

% ROIs for which we extracted peaks of  activation
% Also, those to consider for the expansion intersection
opt.roiList = {'VWFAfr', 'VWFAbr', 'lLO', 'rLO'};

% Radius of the sphere around the peak
opt.radius = 10; 

% Number of voxels in the case of expanding ROI
opt.numVoxels = 150;

% Number of voxels in the case of expanding ROI
opt.numLanguageVoxels = 80;

% Save the ROI?
opt.saveROI = true;

% Specify space accordingly to source of data 
% - IXI549Space for bidspm preprocessing
% - MNI152NLin2009cAsym for fmriprep
% - individual
opt.space = 'MNI152NLin2009cAsym'; 

% description to add to folder name, to distinguish from GLM 
% Possibly obsolete
opt.desc = 'ROI';

% task to analyze
opt.taskName = 'visualLocalizer';


%% Check that everything is in order

% Sphere - all parameters should be present

% Intersection - all parameters should be present

% Expansion - all parameters should be present

% Atlas - all parameters should be present and the atlas should be present


end

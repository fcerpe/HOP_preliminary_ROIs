function opt = roi_option()
% returns a structure that contains the options chosen by the user to run
% createROI (extraction)

if nargin < 1
    opt = [];
end

% suject to run in each group
opt.subjects = {'006','007','008','009','010','011','012','013','018','019','020','021','022','023','024','026','027','028'}; 
% Participants: 
% '006','007','008','009','010','011','012','013','018','019','020','021','022','023','024','026','027','028'

% ROIs for which we extracted peaks of  activation
% Also, those to consider for the expansion intersection
opt.roiList = {'VWFAfr', 'VWFAbr', 'lLO', 'rLO'};

% Radius of the sphere around the peak
opt.radius = 10; 

% Number of voxels in the case of expanding ROI
opt.numVoxels = 115;

% Number of voxels in the case of expanding ROI
opt.numLanguageVoxels = 80;

% Option(s) to execute in vwfa split
% Possible options:
% - 'atlas': VWFA ROI for each subject intersected with aVWFA/pVWFA
%            coordinates form visfatlas, with perVWFA adn lexVWFA from 
%            Lerma-Usabiaga et al. (2018)
% - individual: split each VWFA ROI in anterior / posterior halves
opt.split = {'individual'};

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

% I like chatty outputs
opt.verbosity = 2;

% task to analyze
opt.taskName = 'visualLocalizer';

% PATHS
% The directory where the data are located
opt.dir.root = fullfile(fileparts(mfilename('fullpath')), '..', '..');

opt.dir.raw = fullfile(opt.dir.root, 'inputs', 'raw');
opt.dir.derivatives = fullfile(opt.dir.root, 'outputs', 'derivatives');
opt.dir.preproc = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-preproc');
opt.dir.input = opt.dir.preproc;
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'cpp_spm-rois');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-stats');
opt.dir.cosmo = fullfile(opt.dir.root, 'outputs', 'derivatives', 'CoSMoMVPA');

% Suffix output directory for the saved jobs
opt.jobsDir = fullfile(opt.dir.stats, 'jobs', opt.taskName);
opt.glm.QA.do = false;

% multivariate
opt.pipeline.type = 'roi';

%% DO NOT TOUCH
opt = checkOptions(opt);
saveOptions(opt);

end

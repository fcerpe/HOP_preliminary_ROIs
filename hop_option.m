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

% Avoid possible overlaps for multiple runs, start from scratch everytime
opt = [];
roiMethods = [];


%% Main parameters
% these parameters are valid for all the methods and are stored in 'opt'


% PATHS 
% - root: main folder of the project
% - raw: where your un-preprocessed data is
% - derivatives: main folder where you store the results
% - preproc / stats / roi / cosmo: main folders (subfolders of derivatives)
%                                  to store different steps of the pipeline
% - jobs: where to store the copy of the script 
opt.dir.root = fullfile(fileparts(mfilename('fullpath')), '..', 'VBE_data');
opt.dir.raw = fullfile(opt.dir.root, 'inputs', 'raw');
opt.dir.derivatives = fullfile(opt.dir.root, 'outputs', 'derivatives');
opt.dir.preproc = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-preproc');
opt.dir.rois = fullfile(opt.dir.root, 'outputs', 'derivatives', 'cpp_spm-rois');
opt.dir.stats = fullfile(opt.dir.root, 'outputs', 'derivatives', 'bidspm-stats');
opt.dir.cosmo = fullfile(opt.dir.root, 'outputs', 'derivatives', 'CoSMoMVPA');

opt.dir.output = fullfile(fileparts(mfilename('fullpath')), 'outputs');
opt.dir.jobs = fullfile(fileparts(mfilename('fullpath')), 'jobs');


% DEFAULT OPTIONS 
% Applied in case no specific parameter is given in a method
% (i.e. if you don't choose specific subejcts, will apply the method on all these subs)

% Space of the data 
% Any space is fine, as long as it's consistent with the space
% indicated in the ROI. Will use this parameter to find and to save ROIs
% Examples: 
% - IXI549Space (bidspm preprocessing)
% - MNI152NLin2009cAsym (fmriprep preprocessing)
% - individual
% - MNI
opt.space = 'IXI549Space'; 


% Stats parameters
% used to find correct GLM analysis and necessary files 
% * Task: to find contrasts in the 'intersection' and 'expansion' methods
opt.task = 'visualLocalizer';

% * Node
opt.node = 'localizerGLM';

% * FWHM
opt.fwhm = 6;


% Level on which to draw ROIs
% choose the default option to draw ROIs: 
% - individual: on each subject
% - group
opt.level.default = 'individual';

% In the case of individual ROIs, enter a default list of subjects. If not
% specificed for a give method, will apply the method on subjects from this
% list
opt.level.subjects = {'006','007'}; 


% Save ROIs
% choose whether to save the ROIs and whether to use subfolders for each method
opt.saveROI = true;

% subfolder is only applicable to 'intersection' and 'expansion' methods
% methods 'sphere' and 'atlas' do not include subjects, so they will be
% added to subfolders that specify method used and radius/atlas
opt.saveInSubfolder = false; 


%% Template methods 
% Feel free to modify / copy / delete the following methods to your needs,
% adding all the parameters needed


%% Method #1 - Sphere of [radius]mm in [area]
roiMethods(1).method = 'sphere';

% what radii should the sphere(s) have?
roiMethods(1).radii = [10, 8, 6];

% which areas will we extract? 
roiMethods(1).area = {'VWFA'};

% from which coordinates? 
% (if lateralized, omit controlateral coords using NaNs)
roiMethods(1).coordsL = [-46 -56 -16];
roiMethods(1).coordsR = [NaN NaN NaN];

% for this method, we need a reference image (any beta.nii from any
% subject's GLM). Which path should we follow? 
roiMethods(1).referencePath = fullfile(opt.dir.stats, ['sub-', roiMethods(1).subjects{1}], ...
                                       ['task-', opt.task,'_space-',opt.space,'_FWHM-', num2str(opt.fwhm),'_node-',opt.node], ...
                                       'beta_0001.nii');

% Do you want the ROIs to be merged across hemispheres? 
% Only available in the case you provide two valid coordinates sets for
% left and right hemisphere. If not, won't merge regardless your choice
roiMethods(1).mergeRois = 0;

%% Method #2 - Intersection between [contrast] and [roi]
roiMethods(2).method = 'intersection';

% on which subjects should we extract ROIs?
roiMethods(2).subjects = {'006','007'}; 

% which area will be extracted?
roiMethods(2).area = 'area1';

% details of the ROI
roiMethods(2).roiPath = '../whatever';

% Optional: t-map can be difend either by providing information to extract
% the contrast, or by the exported .nii map
roiMethods(2).tmapPath = '../whatever';
roiMethods(2).task = 'task';
roiMethods(2).contrast = 'contrast';

% which min. number of voxel can we accept? 
roiMethods(2).nbVoxel = 25;


%% Method #3 - Atlas: extraction of [area] from [atlas]
% Assumes that the mask is correct, will trust you on this
roiMethods(3).method = 'atlas';

% atlas and mask names are used to save the mask 
roiMethods(3).atlas = 'atlas';
roiMethods(3).mask = 'label';
roiMethods(3).maskPath = '../whatever';


%% Method #4 - Expansion around [area] in [contrast]
roiMethods(4).method = 'expansion';

% How many voxels we want to reach?
roiMethods(4).nbVoxel = 100;

% Sphere
roiMethods(4).area = {'area1'};
roiMethods(4).coordsL = [0 0 0];
roiMethods(4).coordsR = [0 0 0];

% T-map
roiMethods(4).tmapPath = '../whatever';
roiMethods(4).task = 'task';
roiMethods(4).contrast = 'contrast';



%% Preliminary checks
% is everything in order? 

% Throw an error if something major is missing (paths, areas)
% Add defaults if something minor is missing (subs, radii, voxels)

% Sphere - all parameters should be present
% Intersection - all parameters should be present
% Expansion - all parameters should be present
% Atlas - all parameters should be present and the atlas should be present

end


%% Support functions



















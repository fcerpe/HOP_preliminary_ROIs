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


% Checks on the methods

opt.checks.availableAtlases = {'Brodmann', 'custom'};


%% Template methods 
% Feel free to modify / copy / delete the following methods to your needs,
% adding all the parameters needed


%% Method #1 - Sphere of [radius]mm in [area]
roiMethods(1).method = 'sphere';

% what radii should the sphere(s) have?
% will do all the radii for all the ROI specified
roiMethods(1).radii = [10, 8, 6];

% for this method, we need a reference image (any beta.nii from any
% subject's GLM). Which path should we follow? 
roiMethods(1).referencePath = fullfile(opt.dir.stats, ['sub-', opt.level.subjects{1}], ...
                                       ['task-', opt.task,'_space-',opt.space,'_FWHM-', num2str(opt.fwhm),'_node-',opt.node], ...
                                       'beta_0001.nii');

% Method-specific details for each ROI that you want to create
methodDetails = [];

% - which areas will we extract
% - from which coordinates (if lateralized, leave controlateral empty)
% - on which reference image will we resize the ROI created
% - whether to merge the ROIs or to keep them separate
methodDetails(1).area = 'VWFA';
methodDetails(1).coordsL = [-46 -56 -16];
methodDetails(1).coordsR = [];
methodDetails(1).mergeRois = 0;


% Add all the rois to a struct
roiMethods(1).roisToCreate = methodDetails;


%% Method #2 - Atlas: extraction of [area] from [atlas]
% avaialble atlases are:
% - Broadmann
% - visfatlas (coming soon)
% - fedorenko (coming soon)

roiMethods(2).method = 'atlas';

% for this method, we need a reference image (any beta.nii from any
% subject's GLM). Which path should we follow? 
roiMethods(2).referencePath = fullfile(opt.dir.stats, ['sub-', opt.level.subjects{1}], ...
                                       ['task-', opt.task,'_space-',opt.space,'_FWHM-', num2str(opt.fwhm),'_node-',opt.node], ...
                                       'beta_0001.nii');

% Method-specific details for each ROI that you want to create
methodDetails = [];

% - from which atlas we will extract the ROIs
% - which areas will we extract
% - from which parcels
methodDetails(1).atlas = 'Brodmann';
methodDetails(1).area = 'DLPFC';
methodDetails(1).parcels = [9 46];

methodDetails(2).atlas = 'Brodmann';
methodDetails(2).area = 'LVC';
methodDetails(2).parcels = [17 18];

% are you working on personal / specific parcels? 
% - you can set the atlas to 'custom' and skip parcels
% - you can specify the path to your mask 
methodDetails(3).atlas = 'custom';
methodDetails(3).area = 'AREA';
methodDetails(3).roiPath = 'masks/path-to-custom-parcel';


% Add all the rois to a struct
roiMethods(2).roisToCreate = methodDetails;


%% Method #3 - Intersection between [contrast] and [roi]
roiMethods(3).method = 'intersection';

% on which subjects should we extract ROIs?
% It can be a wildcard 'all', indicating that all subjects in the dataset 
% should be processed
roiMethods(3).subjects = {'006','007'}; 

% which min. number of voxel can we accept? 
roiMethods(3).targetNbVoxels = 25;


% Method-specific details for each ROI that you want to create
methodDetails = [];

% FFA - Faces > Objects
methodDetails(1).area = 'FFA';
methodDetails(1).roiPath = fullfile(opt.dir.output,'radius_10/ROI-FFA_radius-10_space-MNI_resampled-to-sub_binary.nii');
methodDetails(1).tmapPath = [];
methodDetails(1).task = 'loc1';
methodDetails(1).contrast = 'Faces > Objects';

% LOC - Objects > Scrambled
methodDetails(2).area = 'LOC';
methodDetails(2).roiPath = fullfile(opt.dir.output,'radius_10/ROI-LOC_radius-10_space-MNI_resampled-to-sub_binary.nii');
methodDetails(2).tmapPath = [];
methodDetails(2).task = 'loc1';
methodDetails(2).contrast = 'Objects > Scrambled';

% PPA - Scenes > Objects
methodDetails(3).area = 'PPA';
methodDetails(3).roiPath = fullfile(opt.dir.output,'radius_10/ROI-PPA_radius-10_space-MNI_resampled-to-sub_binary.nii');
methodDetails(3).tmapPath = [];
methodDetails(3).task = 'loc1';
methodDetails(3).contrast = 'Scenes > Objects';


% Add all the rois to a struct
roiMethods(3).roisToCreate = methodDetails;


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

% Sphere
% - coords should be empty or triplets

% Intersection

% Expansion

% Atlas
% - atlas specified should be in a list or be 'custom'
% if ~ismember(inputAtlas, supportedAtlases)
%     error('MyComponent:incorrectType',...
%           ['Atlas specified is not supported at the moment. Check the options\n' ...
%           'Supported atlases: Brodmann, custom'])
% end

end


%% Support functions



















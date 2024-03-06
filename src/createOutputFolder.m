function outputFolder = createOutputFolder(opt, m, detail, subject)

    % 04/03/2024 
    % Only works for atlas

    % Creates an output directory for a given subject and ROI if it doesn't already exist.
    %
    % Parameters:
    %   outRoot: The root directory under which all outputs are saved.
    %   outFolder: The folder name specific to the ROI or processing step.
    %   subName: The name of the subject being processed.
    %
    % Example usage:
    %   createOutputDir('/path/to/output/root', 'roi_specific_folder', 'subject_name');
    %
    % This function constructs the output directory path, checks if it exists,
    % and creates it if not, with logging at each step for clarity.

    % if subject argument is passed (i.e. 4 input arguments are passed to the function),
    % make a note to create a subject-specific folder.
    % Otherwise, do not add the subject
    addSub = nargin == 4;

    % get method specifics, if they are not passed as input already
    switch m.method
        case 'sphere'
            methodDetails = ['_radius-', num2str(detail)];
        case 'atlas'
            methodDetails = ['_atlas-', m.atlas];
        case 'intersection'
            methodDetails = ['_joined-', extractIntersection(m)];
        case 'expansion' 
            methodDetails = [];
    end

    % Construct full path
    % - specified (general) output folder
    % - method name and its details
    % - subject if needed (coming soon)
    if addSub, outputFolder = fullfile(opt.dir.output, ['method-', m.method, methodDetails], subject);
    else, outputFolder = fullfile(opt.dir.output, ['method-', m.method, methodDetails]);
    end
    
    % Create the folder, if it does not already exists
    if ~exist(outputFolder, 'dir')

        fprintf('Creating output folder\n');

        % Attempt to create the directory
        [success, msg, msgID] = mkdir(outputFolder);

        if success
            fprintf('Output folder created successfully: %s\n\n', outputFolder);
        else
            error('Failed to create output folder: %s\nMessage ID: %s\n%s\n\n', outputFolder, msgID, msg);
        end

    else
        fprintf('Output folder already exists: ''%s'' \nNo action needed\n\n', outputFolder);
    end

end


%% Subfunctions

% From the path of the ROI used in the intersection between ROI and GLM,
% extract details to summarize the ROI used
% e.g. path = '/method-sphere_radius-10/' -> '10mm-sphere'
%      path = '/method-atlas_atlas-Brodmann/..._ROI-DLPFC' -> 'Brodmann-DLPFC
% other options TBD
function outStr = extractIntersection(m)
    
    % isolate roi name 
    roiParts = split(m.roiPath, '/');
    roiName = roiParts{end}; 

    % split the roi name into different elements
    roiDetails = split(roiName, {'-','_'});

    % find 
    % - 'method' 
    % - details about method (radius, area)
    % - name of the area (label)
    method = roiDetails{find(strcmp('method',roiDetails))+1};
    label = roiDetails{find(strcmp('label',roiDetails))+1};
    switch method
        case 'sphere', detail = roiDetails{find(strcmp('radius',roiDetails))+1};
        case 'atlas', detail = roiDetails{find(strcmp('atlas',roiDetails))+1};
    end

    % compose new name
    switch method
        case 'sphere', folderInfo = [detail, 'mm-sphere-', label];
        case 'atlas', folderInfo = [detail, '-', label];
    end   

    outStr = folderInfo;
    
end



function [tConImg, tStruct] = loadTContrastImage(model, contrastName, GLMdir)
    % Loads T contrast image voxel data for a given contrast name and GLM directory.
    %
    % Parameters:
    %   model: The loaded SPM model object for the current subject and task.
    %   contrastName: The name of the contrast to load.
    %   GLMdir: Directory containing the GLM results for the subject.
    %
    % Returns:
    %   tConImg: The voxel data of the T contrast image.
    %   tStruct: The structure of the T contrast image.
    %
    % Example usage:
    %   [tConImg, tStruct] = loadTContrastImage(model, 'contrastName', '/path/to/GLMdir');
    %
    % This function first attempts to find the specified contrast within the model.
    % If not found, it appends ' - All Sessions' to the name and retries.
    % It then checks for the existence of the T contrast image file and loads its voxel data.

    try
        fprintf('Attempting to load T contrast image for contrast: %s\n', contrastName);
        
        % Attempt to find the specified contrast by name
        t_con = get_contrast_by_name(model, contrastName);
        if isempty(t_con)
            % If not found, try appending ' - All Sessions' and search again
            fprintf('Contrast %s not found, trying with " - All Sessions" suffix.\n', contrastName);
            contrastName = strcat(contrastName, ' - All Sessions');
            t_con = get_contrast_by_name(model, contrastName);
            if isempty(t_con)
                error('ContrastNotFound', 'Cannot find the contrast %s in the design; has it been estimated?', contrastName);
            end
        end

        % Construct the full path to the T contrast image file
        tConFname = fullfile(GLMdir, t_con.Vspm.fname);

        % Verify the existence of the T contrast image file
        if ~exist(tConFname, 'file')
            error('FileNotFound', 'Cannot find T image %s; has it been estimated?', tConFname);
        else
            fprintf('T contrast image found: %s\n', tConFname);
        end

        % Load the voxel data from the T contrast image
        fprintf('Loading voxel data from T contrast image...\n');
        tStruct = spm_vol(tConFname);
        tConImg = spm_read_vols(tStruct);
        fprintf('Voxel data loaded successfully.\n');

    catch ME
        switch ME.identifier
            case 'ContrastNotFound'
                fprintf('Error: %s\n', ME.message);
            case 'FileNotFound'
                fprintf('Error: %s\n', ME.message);
            otherwise
                fprintf('An unexpected error occurred: %s\n', ME.message);
        end
        tConImg = []; % Return empty in case of error
        tConFname = '';
    end
end
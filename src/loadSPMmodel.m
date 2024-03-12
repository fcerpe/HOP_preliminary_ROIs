function [model, loadSuccess] = loadSPMmodel(GLMdir, subName, taskName)
    % Attempts to load the SPM model for a given subject and task.
    %
    % Parameters:
    %   GLMdir: The directory containing the GLM results for the subject.
    %   subName: The name of the current subject being processed.
    %   taskName: The name of the task for which the model is being loaded.
    %
    % Returns:
    %   model: The loaded SPM model object, or empty if loading failed.
    %   loadSuccess: A boolean indicating whether the model was successfully loaded.
    %
    % Example usage:
    %   [model, loadSuccess] = loadSPMModel('/path/to/GLMdir', 'subject1', 'task1');
    %
    % This function tries to load the SPM.mat file and handles any errors that occur,
    % logging appropriate messages and returning a status indicator.

    fprintf('Attempting to load SPM model for subject %s, task %s...\n', subName, taskName);
    model = []; % Initialize model as empty
    loadSuccess = false; % Initialize success status as false

    try
        % Construct the path to the SPM.mat file
        spmPath = fullfile(GLMdir, 'SPM.mat');

        % Attempt to load the SPM model
        model = mardo(spmPath);
        
        % If successful, set the success status to true
        loadSuccess = true;
        fprintf('SPM model loaded successfully.\n');

    catch ME
        % Handle errors that occur during model loading
        fprintf('WARNING: Error loading SPM model for %s, task %s: %s\n', subName, taskName, ME.message);
        fprintf('Skipping to the next iteration.\n');
        % No need to set model or loadSuccess as they are already initialized to their failure states
    end

    % Return the model (empty if failed) and the success status
    return
end
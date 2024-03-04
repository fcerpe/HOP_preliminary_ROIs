% Import BA atlas and gunzip if needed

% If a valid atlas is provided, will extract (unzip, load, get volumes) it.
% If a custom parcels is provided, load (maybe unzip) it as it is
% atlas = importAtlas(m)
function atlasPath = importAtlas(m)

% Import specified atlas 
% 
% If present in the list, import it
% otherwise take whatever is given as input


supportedAtlases = {'Brodmann', 'custom'};


% Check that the atlas requested is present among the atlases implemented
if ismember(m.atlas, supportedAtlases) 

    % extract atlas path
    switch m.atlas 
        case 'Brodmann'
            atlasPath = 'masks/Brodmann/brodmann.nii.gz';

        % Treat the custom atlas as a regular one, but specify the path to
        % be the one provided as parameter
        case 'custom'
            atlasPath = m.roiPath;
    end

% The parameter passed as atlas is not implemented yet
else
    error('MyComponent:incorrectType',...
          ['Atlas specified is not supported at the moment. Check the options\n' ...
          'Supported atlases: Brodmann, custom'])

end


end
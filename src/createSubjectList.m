function list = createSubjectList(opt, m)

% Create a subject list to pick the requested subjects
% Not perfect, needs control over what is passed as m.subjects

% Check whether method's subjects are empty, in case use default ones
if isempty(m.subjects)
    list = opt.level.subjects;

% Check whether method's subjects are 'all', in case use every 'sub-*'
% present in the stats / GLM folder
elseif strcmp(m.subjects,'all')
    list = opt.level.subjects;

% Check whether method's subjects are specified, in case use those
else
    list = m.subjects;

end


end
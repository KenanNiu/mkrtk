function Loa = tendon_calculator(handles,studyname)
% A work in progress to calculate tendon lines of action
%
% This function handles interaction with the gui data to choose which
% method of calculation to use, call that method, and do the associated
% plotting

% The process for using this would be:
%   1) Load models in MRIMagic
%   2) Register 3D to dynamic
%   3) Smooth motion
%   4) Define / calculate helical axes
%   5) Run this function...

% This section is mere convenience for debugging:
if nargin == 0
    hf = findall(0,'Name','MRIMagic');
    handles = guidata(hf);
end

fprintf(2,'1. You will probably want to make sure your motion is smoothed before running this function.\n')
fprintf(2,'2. Make sure you have a helical axis defined\n');
fprintf(2,'3. Define a line of action through the normal GUI (which doesn''t yet run calcs)\n\n');
fprintf(2,'4. Finally, call this function\n\n');

assert(~isempty(handles.LineOfAction),'You must define a line of action')

% Get model of the tendon which defines the line of action:
tendon = handles.Models(strcmpi({handles.Models.Tag},handles.LineOfAction(1).Item));
%tendon = get_tendon_model(handles.Models,handles.LineOfAction)


joint = joint_name({handles.Models.Tag});

switch joint
    case 'ankle'
        Loa = tendon_loa(tendon,studyname,'ankle');
        disp(' ')
        disp(studyname)
        
    case 'knee'
        Loa = tendon_loa(tendon,studyname,'knee');
        
end
        
%[~,home] = system('echo $HOME'); home = deblank(home);
%print(hf,'-dtiff', [home filesep 'Desktop' filesep studyname '.tiff'])


% ------------------------------------------------------------------------
function joint = joint_name(tags)
% Get the joint name.  Return either:
%   'ankle'
%   'knee'
%
% Input TAGS comes from handles like this:
% tags = {handles.Models.Tag}

tags = lower(tags);

ankle_keys = {'tal','cal','ach'};
knee_keys  = {'fem'};

for j = 1:numel(ankle_keys)
    tf = ~cellfun(@isempty,strfind(tags,ankle_keys{j}));
    if any(tf)
        joint = 'ankle';
        return
    end
end

for j = 1:numel(knee_keys)
    tf = ~cellfun(@isempty,strfind(tags,knee_keys{j}));
    if any(tf)
        joint = 'knee';
        return
    end
end
warning('Could not determine which joint we''re looking at')
disp('**** -- Selecting Ankle for default -- ****')
joint = 'ankle';




% ------------------------------------------------------------------------
function mdl = get_tendon_model(Models)
names = {Models.Tag};
lookfor = @(tag)strfind(lower(names),tag);
locs = lookfor('ach');

if all(cellfun(@isempty,locs))
    locs = lookfor('ten');
end
if all(cellfun(@isempty,locs))
    error('Could not find which model is the achilles tendon.  Rename model or troubleshoot')
end

model_id = ~cellfun(@isempty,locs);
mdl = Models(model_id);


function Loa = tendon_calculator(handles)
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

fprintf(2,'Make sure your motion is smoothed before defining the helical axis!!!\n');
fprintf(2,'1. Change smoothing parameters\n');
fprintf(2,'2. Smooth the motion: Dev>Smooth\n');
fprintf(2,'3. Define / re-define / calculate helical axis\n');
fprintf(2,'4. Finally, call this function\n\n');

assert(~isempty(handles.LineOfAction),'You must define a line of action')

% Get model of the tendon which defines the line of action:
tendon = handles.Models(strcmpi({handles.Models.Tag},handles.LineOfAction(1).Item));
%tendon = get_tendon_model(handles.Models,handles.LineOfAction)

studyname = study_name_from_path(tendon.LoRes(1).Path);

joint = joint_name({handles.Models.Tag});

switch joint
    case 'ankle'
        Loa = tendon_loa(tendon,studyname,'ankle');
        disp(' ')
        disp(studyname)
        ankle_plotting(handles,Loa)
        
    case 'knee'
        Loa = tendon_loa(tendon,studyname,'knee');
        %knee_plotting(handles,Loa)
        ankle_plotting(handles,Loa)
        
end
        
[~,home] = system('echo $HOME'); home = deblank(home);
%print(hf,'-dtiff', [home filesep 'Desktop' filesep studyname '.tiff'])


% ------------------------------------------------------------------------
function ankle_plotting(handles,Loa)

Hax = handles.HelicalAxis(end);
theta = handles.HelicalAxis(end).Angle*180/pi;
Hax = [Hax.Point Hax.Axis];

n = size(Loa,1);
Marm = NaN(n,1); 

% Now the two data sets are in the same form:
for j = 1:n
    p1 = Loa(j,1:3);
    p2 = p1 + Loa(j,4:6);
    p3 = Hax(j,1:3);
    p4 = p3 + Hax(j,4:6);
    
    Marm(j) = momentArm(p1,p2,p3,p4);
end

studyname = study_name_from_path(handles.Models(1).LoRes(1).Path);
studyname = regexprep(studyname,{' ','*','(',')','\.','&','\$','@','#','\\','\/','\^'},'_');
var_name = [studyname '_Moment_Arm'];
assignin('base',var_name,Marm)
evalin('base',var_name)
clipboard('copy',sprintf('%10.4f\n',Marm))
disp('|--> Copied to clipboard')

%
hf = figure;
ax1 = subplot(2,1,1);
plot(Marm,'b-')
hold on,
plot(Marm,'r*')
%set(ax1,'Ycolor',[0 0 0.8])
xl = xlim(ax1);
yl = ylim(ax1);
xlim(ax1,[0 xl(2)])
ylim(ax1,[0 ceil(yl(2)*5)/5])
xlabel(ax1,'Phase ID')
ylabel(ax1,'Moment arm [mm]')
grid on

ax2 = subplot(2,1,2);
plot(ax2,theta,'Color',[0 0.9 0]), hold on
plot(ax2,zeros(size(theta)),'color',[0 0.5 0],'LineStyle','--')
%set(ax2,'Color','none')
%set(ax2,'YColor',[0 0.5 0])
%set(ax2,'YAxisLocation','right')
%n = numel(get(ax1,'YTick'));
%tlim = round(max(abs(theta)*n));
%set(ax2,'Ylim',[-tlim n*tlim])
%set(ax2,'YTick',linspace(-tlim,tlim,n))
ylim(ax2, ylim(ax2)*2)
ylabel(ax2,'Helical axis rotation angle [deg]')
grid on


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
error('Could not determine which joint we''re looking at')



% ------------------------------------------------------------------------
function name = study_name_from_path(str)
% Get the study name for the pilot studies:
snames = {'pilot', 'mgr'};
for j = 1:numel(snames)
    snidx = strfind(lower(str),snames{j});
    if isempty(snidx)
        continue
    end
    fsidx = strfind(lower(str),filesep);
    if isempty(fsidx)                           %
        if filesep == '/'                       % Handle strings from another
            fsep = '\';                         %  platform. 
        else                                    %
            fsep = '/';                         %
        end                                     %
        fsidx = strfind(lower(str),fsep);       %
    end                                         %
    stridx = fsidx([find(fsidx < snidx(1),1,'last'), find(fsidx > snidx(1),1,'first')]);
    name = str(stridx(1)+1:stridx(2)-1);
    if ~isempty(name)
    end
end

assert(~isempty(name))


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


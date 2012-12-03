function varargout = Registration(varargin)
% REGISTRATION MATLAB code for Registration.fig
%      REGISTRATION, by itself, creates a new REGISTRATION or raises the existing
%      singleton*.
%
%      H = REGISTRATION returns the handle to a new REGISTRATION or the handle to
%      the existing singleton*.
%
%      REGISTRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTRATION.M with the given input arguments.
%
%      REGISTRATION('Property','Value',...) creates a new REGISTRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Registration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Registration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Registration

% Last Modified by GUIDE v2.5 22-Nov-2012 15:22:24

% ------------------------------------------------------------------------
% NOTES 
%
%  - Gaussian Mixture Model:
%    http://www.mathworks.com/matlabcentral/fileexchange/20856-robust-point-set-registration-using-mixture-of-gaussians    
%
%  - Build notes for LIBICP: http://www.cvlibs.net/software/libicp.html
%       - Download LIBICP from above link & extract
%       - Download BOOST libraries from www.boost.org & extract
%       - Place the folder /boost_1_49_0/boost/ into the /libicp/src/ directory
%       - Edit the first line of icpMex.cpp to use <> instead of "", as
%           follows:
%               #include <mex.h>
%       - Run make.m from the /libicp/matlab/ directory
%
%  - Point Cloud Library (PCL):
%       - http://pointclouds.org/
%
%  - CompareCloud:
%       - http://www.danielgm.net/cc/
%
% NEW FEATURES TO ADD:
%
%   - "Settings" interface
%       - Persistent sessions?  Or just remember paths
%       - Enable/disable parallel processing with user selection of # cores
%           - # cores Mac: 
%               [retval,ncores] = system('sysctl -n hw.ncpu')
%           - # cores Windows:
%               ?[retval,ncores] = system('echo %NUMBER_OF_PROCESSORS%')
%   
% CHANGES / IMPROVEMENTS:
%
%   - Re-design AnalysisGui
%       - Stability: on load, ensure reference items exist, if not, drop it
%   - ModelLoader
%       - Add Drag & Drop support for adding *.mat files (see uitree.m)
%   - Registration
%       - Include outlier rejection according to Masuda '96
%       - Use LMS estimator (see Masuda '96) instead of RMS error
%   - Documentation
%       - Latex user manual / quick start guide.
%
% KNOWN BUGS:
%   - Updating the slider min/max isn't clear
%       - leftovers from combined program (slider visiblity listener)
%       - The following call probably has no effect now:
%           PhaseSlider_Callback(handles.PhaseSlider, 'init', handles)
%   - Analysis GUI is incomplete / has bugs.
%   
% 
% ------------------------------------------------------------------------

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Registration_OpeningFcn, ...
                   'gui_OutputFcn',  @Registration_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%=========================================================================
%               ESSENTIAL GUI FUNCTIONS
%=========================================================================
% --- Executes just before Registration is made visible.
function Registration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Registration (see VARARGIN)

% Prevent re-launch:
if isfield(handles,'initialized')
    return
else
    handles.initialized = 1;
end

% Mode configuration for helper functions:
modestr = '3d';

% Include some important functions:
configurePaths(modestr)
check_dependencies(modestr);        % Check all external dependencies exist

% Disble work in progress:
if ~isdeveloper
    set(handles.Menu_Dev,'Enable','off')
end

% Choose default command line output for Registration
handles.output = hObject;

% Configure context menus:
handles.hscmenu  = createContextMenu('staticCloud');    % Cloud context menu

% 3D Annotations:
handles.AnnotationCanvas3D = annotationManager(...  %\_ Create a canvas (axes)
    handles.axes1,'NewCanvas','AnnotationCanvas3D');%/
handles.PhaseText = annotationManager(...
    handles.AnnotationCanvas3D,'PhaseText','SouthWest');
set(handles.PhaseText,...
    'Color',[0 0 0],...
    'Units','pixels','Position',[5,20])

% Set up toolbars:
handles = initTools(handles,modestr);

% Initialise our session variables:
handles = Session.new(handles);

% Reserve space for axis definitions:
%handles.HelicalAxis  = [];
%handles.LineOfAction = [];
%handles.MomentArm    = [];

% Configure view mode:
configureView3(handles)

% Position gui nicely:
movegui(handles.figure1,'north')

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = Registration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Slider width/position:
spsn = get(handles.PhaseSlider,'Position');
fpsn = get(handles.figure1,'Position');
set(handles.PhaseSlider,'Position',[0 0 fpsn(3)-15 spsn(4)])


%=========================================================================
%               GUI USER FUNCTIONS
%=========================================================================


% ------------------------------------------------------------------------
function Menu_File_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ------------------------------------------------------------------------
function MI_OpenSession_Callback(hObject, eventdata, handles)
% hObject    handle to MI_OpenSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Lock figure:
locker = FigLocker.Lock(handles.figure1);

% Get file from usesr
[filename, pathname] = uigetfile(...
    {'*.mat','MAT-files (*.mat)'},'Load session',handles.sessionPath);

if isequal(filename,0)
    locker.unlock;
    return
end

locker = FigLocker.Lock(handles.figure1);

[handles,ok] = Session.load(handles,[pathname filename]);
if ok
    % Update handles:
    guidata(hObject,handles)
    
    % Configure/reset the view:
    configureView3(handles)
end

locker.unlock;


% ------------------------------------------------------------------------
function MI_SaveSession_Callback(hObject, eventdata, handles)
% hObject    handle to MI_SaveSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Lock figure
locker = FigLocker.Lock(handles.figure1);

% Default path:
seedPath = handles.sessionPath;
if isempty(seedPath) && ~isempty(handles.userPath)
    seedPath = handles.userPath;
end
% Shared function:
[pathname,ok] = Session.save(seedPath,handles);
if ok
    handles.sessionPath = pathname;
    guidata(hObject,handles)
end

% Unlock figure:
locker.unlock;


% ------------------------------------------------------------------------
function MI_NewSession_Callback(hObject, eventdata, handles)
% hObject    handle to MI_NewSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = Session.reset('soft',handles);
if ~isempty(handles)
    configureView3(handles)
end

% ------------------------------------------------------------------------
function MI_Quit_Callback(hObject, eventdata, handles)
% hObject    handle to MI_Quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)


% ------------------------------------------------------------------------
function Menu_Dev_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Dev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ------------------------------------------------------------------------
function Menu_Actions_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Actions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ------------------------------------------------------------------------
function MI_LoadDisplayModels_Callback(hObject, eventdata, handles)
% hObject    handle to MI_LoadDisplayModels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Lock the figure while we're busy with the other GUI:
locker = FigLocker.Lock(handles.figure1);

% Launch the GUI, with error protection:
try
    hf = ModelLoader(handles.Models);
    % Add listener - When window closes, HitTest gets turned off, which
    % triggers this listener, which does all the data & graphics updates for
    % the loading of models.
    addlistener(hf,'HitTest','PreSet',@modelLoadUpdater);
    
catch ME
    % Something went wrong...
    locker.unlock;
    rethrow(ME)
end



% --------------------------------------------------------------------
function MI_Registration_Callback(hObject, eventdata, handles)
% hObject    handle to MI_Registration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

launch_Solver_gui(handles)



% --------------------------------------------------------------------
function MI_Kinematics_Callback(hObject, eventdata, handles)
% hObject    handle to MI_Kinematics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

launch_Analysis_gui(handles)


% ------------------------------------------------------------------------
function PhaseSlider_Callback(hObject, eventdata, handles)
% hObject    handle to PhaseSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Set up slider stops:
nmin = get(hObject,'Min');
nmax = get(hObject,'Max');

% Limit the value to an index in range:
limit = @(val,vmin,vmax)max([vmin min([vmax val])]);
ind = round(get(hObject,'Value'));
ind = limit(ind,nmin,nmax);

% Modify value only if it has been changed in previous lines
if ~isequal(ind,get(hObject,'Value'))
    
    set(hObject,'Value',ind)
end

% Update Text:
annotationManager(handles.PhaseText,ind,nmax)

% Update the phase-based display:
phaseDisplay(handles);


% ------------------------------------------------------------------------
function PhaseSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PhaseSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% Configure / check the slider every time it becomes visible:
addlistener(hObject,'Visible','PreSet',@configurePhaseSlider);


% ------------------------------------------------------------------------
function MI_Replay_Callback(hObject, eventdata, handles)
% Replay motion through all phases.
% See also MI_ExportMovie_Callback, which does a similar thing to create a
% video.

dt = 0.2;  % Target frame time (1/fps)

np = get(handles.PhaseSlider,'Max');

% Get axis limits for the whole cycle:
[xlims,ylims,zlims] = get_bounding_axis_limits(handles);

for j = 1:np
    t = tic;
    % Update display:
    set(handles.PhaseSlider,'Value',j);     % Set phase
    runCallback(handles.PhaseSlider);       % Force draw
    % Set axis limits
    set(handles.axes1,'Xlim',xlims)
    set(handles.axes1,'Ylim',ylims)
    set(handles.axes1,'Zlim',zlims)
    % Draw & pause
    drawnow
    et = toc(t);    
    pause(dt-et);
end


% ------------------------------------------------------------------------
function MI_ExportMovie_Callback(hObject, eventdata, handles)
% hObject    handle to MI_ExportMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % Configure video file:
% [fname,pname] = uiputfile({'*.avi','AVI File (*.avi)'},...
%     'Save video file as',[handles.userPath 'movie.avi']);     
% vidObj = VideoWriter([pname fname],'Motion JPEG AVI');
% vidObj.FrameRate = 20;
% vidObj.Quality = 100;
% open(vidObj);
% set(handles.axes1,'Visible','off')
% set(handles.axes1,'CameraViewAngleMode','manual')
% delta = 2;
% for alfa = 0:delta:360-delta
%     camorbit(handles.axes1,delta,0)
%     drawnow
%     currFrame = getframe(handles.figure1);
%     writeVideo(vidObj,currFrame);
% end
% close(vidObj);
% set(handles.axes1,'Visible','on')
% set(handles.axes1,'CameraViewAngleMode','auto')
% return

% Make a movie of the motion through all phases

% Get axis limits for the whole cycle:
[xlims,ylims,zlims] = get_bounding_axis_limits(handles);
    
% Configure video file:
[fname,pname] = uiputfile({'*.avi','AVI File (*.avi)'},...
    'Save video file as',[handles.userPath 'movie.avi']);     
vidObj = VideoWriter([pname fname],'Motion JPEG AVI');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);

% Record frames:
np = get(handles.PhaseSlider,'Max');
for j = 1:np
    % Set the phase:
    set(handles.PhaseSlider,'Value',j)
    runCallback(handles.PhaseSlider)
    
    % Set axis limits
    set(handles.axes1,'Xlim',xlims)
    set(handles.axes1,'Ylim',ylims)
    set(handles.axes1,'Zlim',zlims)
    
    % Get the frame:
    drawnow
    currFrame = getframe(handles.axes1);
    writeVideo(vidObj,currFrame);
end

close(vidObj);


% ------------------------------------------------------------------------
function MI_CopyFigure_Callback(hObject, eventdata, handles)
% hObject    handle to MI_CopyFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hf = figure('Visible','off');
copyobj(handles.axes1,hf)
set(hf,'Visible','on')


% ------------------------------------------------------------------------
function MI_PublishHandles_Callback(hObject, eventdata, handles)
% hObject    handle to MI_PublishHandles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','handles',handles)
evalin('base','handles')



% ------------------------------------------------------------------------
function MI_SmoothMotion_Callback(hObject, eventdata, handles)
% hObject    handle to MI_SmoothMotion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mdls = handles.Models;

if all( cellfun(@isempty,{mdls.q}) )    % Need to have motion data to work on
    disp('Nothing to do...')
    return
end

StateSmoother;


% ------------------------------------------------------------------------
function MI_EnableCursor_Callback(hObject, eventdata, handles)
% hObject    handle to MI_EnableCursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

show_hide_edit_cursor(hObject)



% --------------------------------------------------------------------
function MI_TendonLoA_Callback(hObject, eventdata, handles)
% hObject    handle to MI_TendonLoA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


addpath('tendon_loa_wip')
Loa_data = tendon_calculator(handles);

handles.LineOfAction(1).Point  = Loa_data(:,1:3);
handles.LineOfAction(1).Vector = Loa_data(:,4:6);

guidata(hObject,handles)


% ------------------------------------------------------------------------
function MI_Test_Callback(hObject, eventdata, handles)
% hObject    handle to MI_Test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Nothing here...')

%{
handles.MomentArm = [];
handles.LineOfAction = [];
handles.HelicalAxis = [];
guidata(hObject,handles)
return
%}

%{
% Plot the quaternion components:
for m = 1:numel(handles.Models)
    qd = handles.Models(m).q.unwrap.double;
    figure('Name',handles.Models(m).Tag)
    subplot(2,1,1)
    plot(qd(1,:))
    title('Real component')
    subplot(2,1,2)
    plot(qd(2:end,:)')
    legend('i','j','k')
%     
%     figure
%     ttl = {'Angle','i','j','k'};
%     for j = 1:4
%         ax(j) = subplot(4,1,j);
%         plot(qd(j,:))
%         title(ttl{j})
%     end
end

%}

%


%=========================================================================
%               TOOLBAR CALLBACKS & FUNCTIONS
%=========================================================================



% ------------------------------------------------------------------------
function CloudVisibility_Callback(hObject, eventdata, handles)
% Funtion to show / hide the static / dynamic clouds in positions specified
% by PhaseSlider, after running ICP calculations to determine those positions.

% Refresh display with current toggle states:
phaseDisplay(handles);


% ------------------------------------------------------------------------
function ClearAnalysis_Callback(hObject, eventdata, handles)
% Function to erase handles.motion

if isempty(handles.Models)
    return
end

btn = questdlg('Clicking OK will erase the motion reconstruction. Are you sure?',...
    'Erase motion calculations',...
    'Cancel','OK','Cancel');

if ~strcmpi(btn,'OK')
    return
end

% Remove these fields if they exist:
if isfield(handles.Models,'xraw')
    handles.Models = rmfield(handles.Models,'xraw');
end
if isfield(handles.Models,'qraw')
    handles.Models = rmfield(handles.Models,'qraw');
end

% And empty these fields:
[handles.Models.q] = deal([]);
[handles.Models.x] = deal([]);

% Update guidata
guidata(hObject,handles)   

% Refresh display:
phaseDisplay(handles)


% ------------------------------------------------------------------------
%function HelicalAxis_Callback(hObject, eventdata, handles)
%launch_Analysis_gui(handles)


% ------------------------------------------------------------------------
%function ICPCalculateTool_Callback(hObject, eventdata, handles)
%launch_Solver_gui(handles)

% ------------------------------------------------------------------------
function ShowAxesTool_Callback(hObject, eventdata, handles)

% Refresh display with current toggle states:
phaseDisplay(handles);

%=========================================================================
%               CUSTOM HELPER FUNCTIONS
%=========================================================================

% ------------------------------------------------------------------------
function configureView3(handles)
% CONFIGUREVIEW3 Configure the view.
%
% This function is called when a new set of models are loaded, or when a
% new session is loaded
%
% See also PHASEDISPLAY

hf = handles.figure1;   % Shorthand
ha = handles.axes1;

set(hf,'Color',[1 1 1]*0.75)
set(ha,...
    'Color','none',...
    'Units','normalized',...    
    'Position',[0.05 0.1 0.9 0.85])
%axis(handles.axes1,'vis3d')     % Equal aspect ratios
%view(handles.axes1,3)           % Default 3D camera view
%grid(handles.axes1,'on')
hold(ha,'on')
grid(ha,'on')
view(ha,3)
axis(ha,'equal','tight')

configurePhaseSlider([],handles.PhaseSlider)  % Update slider limits/steps
figure1_ResizeFcn(hf,[],handles)        % Adjust slider size/position

% Make sure something is being displayed:
if all(strcmpi(get(...
        [handles.ShowStaticTool handles.ShowDynamicTool],'State'),'off'))
    set(handles.ShowStaticTool,'State','on')
end
    
phaseDisplay(handles)

% Axes properties:
%set(ha,'Position',[0.1 0.1 0.8 0.8])
%axis(ha,'equal','tight');   % Make axes tight
%grid(ha,'on');              % Grid on
%set(ha,'Color','none')      % No background colour


% ------------------------------------------------------------------------
function configurePhaseSlider(~,arg2)
% Configure handles.PhaseSlider
% ARG2 could be an event (when called by a listener) or a handle (when
% called directly in code)

% Get handles:
if isa(arg2,'handle.PropertySetEventData')
    handles = guidata(arg2.AffectedObject);
    hObject = handles.PhaseSlider;
else
    hObject = arg2;
    handles = guidata(hObject);
end

if ~isfield(handles,'Models') || isempty(handles.Models)
    set(hObject,'Visible','off')
    return
end

% Set up slider stops & steps:
nmin = 1;
nmax = max( cellfun(@numel,{handles.Models(:).LoRes}) );
if isempty(nmax) || nmax <= nmin
    return
elseif nmax < 10
    incr = [1 2];
else
    incr = [1 5];
end

% Limit the value to an index in range:
limit = @(val,vmin,vmax)max([vmin min([vmax val])]);
ind = round(get(hObject,'Value'));
ind = limit(ind,nmin,nmax);

% Configure slider:
steps = incr/(nmax-nmin);
set(hObject,...
    'Min',nmin,...
    'Max',nmax,...
    'Value',ind,...
    'SliderStep',steps)

% Update phase text:
annotationManager(handles.PhaseText,ind,nmax)


% ------------------------------------------------------------------------
function [xlims,ylims,zlims] = get_bounding_axis_limits(handles)
% Determine the axis limits which bound all objects displayed across all
% phases.
% The easiest way of doing this is to play through the sequence, collect
% axis limits at each phase, then run minx() & max() for each coordinate.

np = get(handles.PhaseSlider,'Max');

Xlims = NaN(np,2);
Ylims = NaN(np,2);
Zlims = NaN(np,2);

% First pass: get suitable axis limits
for j = 1:np
    % Set the phase:
    set(handles.PhaseSlider,'Value',j)
    runCallback(handles.PhaseSlider)
    
    % First pass: get axis limits
    Xlims(j,:) = get(handles.axes1,'Xlim');
    Ylims(j,:) = get(handles.axes1,'Ylim');
    Zlims(j,:) = get(handles.axes1,'Zlim');
end

% Get axis limits:
xlims = [min(Xlims(:,1)) max(Xlims(:,2))];
ylims = [min(Ylims(:,1)) max(Ylims(:,2))];
zlims = [min(Zlims(:,1)) max(Zlims(:,2))];

% ------------------------------------------------------------------------
function launch_Analysis_gui(handles)
% Function to launch the analysis interface

mtags = {handles.Models.Tag};

% For helical axis definitions, use only models that have both static and
% dynamic clouds:
sok = ~cellfun(@isempty,{handles.Models.HiRes});
dok = ~cellfun(@isempty,{handles.Models.LoRes});
ok = sok & dok;
haxCandTags = mtags(ok);

% % Ensure we have some:
% if isempty(mtags)
%     msg = {'There are no models between which helical axes can be defined.';
%         '';
%         ['Make sure you have some models loaded, and that they have both ',...
%         'static and dynamic clouds, then try again.']};
%     warndlg(msg,'Unable to define helical axes','modal')
%     return
% end  

% For line of action definitions, we could use only models that have only
% dynamic clouds, but let's just populate with all for now:
loaCandTags = mtags;

% Launch the GUI, with error protection:
hf = Analysis3DGui(haxCandTags, loaCandTags,...
    handles.HelicalAxis,...
    handles.LineOfAction,...
    handles.MomentArm);
%set(hf,'WindowStyle','modal');  % Optional

% Add listener - When window closes, the HitTest property of ModelLoader
% gets turned off, which triggers this listener which does all the data &
% graphics updates for the loading of models.
addlistener(hf,'HitTest','PreSet',@analysisUpdater);


% ------------------------------------------------------------------------
function launch_Solver_gui(handles)
% Function to launch the registration solver

% Check if we have files:
if isempty(handles.Models)
    warndlg('No files loaded','','modal');
    return
elseif all( cellfun(@isempty,{handles.Models.LoRes}) )
    warndlg('No dynamic files loaded','','modal')
    return
end

if isfield(handles.Models,'qraw')
    for j = 1:numel(handles.Models)
        handles.Models(j).q = handles.Models(j).qraw;
        handles.Models(j).x = handles.Models(j).xraw;
    end
    handles.Models = rmfield(handles.Models,{'qraw','xraw'});
end

% Lock the figure:
locker = FigLocker.Lock(handles.figure1);

% Open the Gui, with error protection:
try
    p = get(handles.PhaseSlider,'Value');      % Current phase
    hf = Solver(handles.Models,p);
    addlistener(hf,'HitTest','PreSet',@solverUpdater);
catch ME
    % Something went wrong... probably input arguments.
    locker.unlock;
    rethrow(ME)
end


% ------------------------------------------------------------------------
function [handles,guiOutput] = getGuiOutput(event)
%GETGUIOUTPUT Get output from GUI which is closing, and handles to main GUI

% Get the handle to the main GUI (Registration)
%   (This function will not receive handles as an input, neither will the
%   main invoking function be in the debug stack because it is invoked by
%   the closing function of the GUI which is closing, so we have to do it
%   another way):  
hmain = findobj('-regexp','Filename',[mfilename '.fig']);   % Get the fig associated with this file
handles = guidata(hmain);                                   % Get it's guidata

% Now get the outputs of the assisting GUI, which has actually invoked this
% function:
loaderhdls = guidata(event.AffectedObject);
guiOutput = loaderhdls.output;


% ------------------------------------------------------------------------
function modelLoadUpdater(~,event) 
%MODELLOADUPDATER Listener callback for getting output from ModelLoader
%
% This function is a listener callback which gets invoked during the window
% closing process of ModelLoader.  The listener needs to be added by any
% function which invokes ModelLoader and this function gets the output of
% that GUI and saves/updates all the info in the handles of Registration (on
% successful close)

% Get handles of main GUI & output from closing GUI
[handles,data] = getGuiOutput(event);

% Unlock the GUI:
FigLocker.Unlock(handles.figure1);

% Check output & store the return data (if we have it):
if isa(data,'Bone')
    
    % We have got a valid output
    handles.Models = data;
    
    % Check button states
    if all(cellfun(@isempty,{handles.Models.LoRes}))
        set(handles.ShowDynamicTool,'State','off')
    end
    if isequal('off',get(handles.ShowStaticTool,'State')) &&...
            isequal('off',get(handles.ShowDynamicTool,'State'))
        % Make sure something is being shown:
        set(handles.ShowStaticTool,'State','on')
    end
    
    % Number of objects:
    no = numel(handles.Models);
    
    % Update the working path:
    if no > 0
        pset = {handles.Models(1).HiRes(:).Path, ... % by concatenating these
            handles.Models(1).LoRes(:).Path};        % we'll ensure we get a path
        pth = pset{1};
        handles.userPath = ...                    % Update user path with
            pth(1:find(pth == filesep,2,'last')); % one up from parent directory
    end
    
    % Update the guidata
    guidata(handles.figure1,handles)
    
    % Check ColorOrder:
    cOrder = get(handles.axes1,'ColorOrder');
    if no > size(cOrder,1)
        newcOrder = lines(no);                    % Create new list, large enough
        newcOrder(1:size(cOrder,1),:) = cOrder;   % Fill with old colours as subset
        set(handles.axes1,'ColorOrder',newcOrder) % Update
    end
    
    % Configure the slider:
    configurePhaseSlider([],handles.PhaseSlider)
    
    % Now we have to update the display:
    phaseDisplay(handles);
    
elseif isempty(data)
    % Cases:
    %   - Figure closed with "X"
    %   - Figure closed with "Cancel"
    %   - Figure closed with "OK" and
    % Do nothing
    return
    
else
    error('Unhandled output from ModelLoader')
    
end


% ------------------------------------------------------------------------
function solverUpdater(~,event)

% Get handles of main GUI & output from closing GUI
[handles,models] = getGuiOutput(event);

hl = FigLocker.GetLocker(handles.figure1);

% If data is empty, user cancelled:
if isempty(models)
    fprintf(2,'Changes discarded...\n');
    try hl.unlock; end %#ok<*TRYNC>
    return  % Do nothing
end

% We will update the models, but first lets see if there have been any
% changes that cascade into other things:

% Have any of the items referenced by helical axes changed?
tf = haxRefsChanged(handles.HelicalAxis,handles.Models,models);

% Update models
handles.Models = models;

% Now perform any calcs that are needed:
if any(tf)
    hl.settext('Updating helical axes...')
    hl.setprogress(inf)
    handles.HelicalAxis(tf) = calcHelicalAxes(handles.HelicalAxis(tf),models);
end

% Now recalculate helical axes if required:

handles.HelicalAxis = calcHelicalAxes(handles.HelicalAxis,handles.Models);

runCallback(handles.PhaseSlider, 'init')
guidata(handles.figure1,handles)

hl.unlock;


% ------------------------------------------------------------------------
function tf = haxRefsChanged(haxes,old_models,new_models)
%HAXREFSCHANGED Test to see if the data that helical axes depend upon has
%been changed.  
tf = false(size(haxes));
mtags = {old_models.Tag};
for j = 1:numel(haxes)
    m1id = strcmp(haxes(j).Item1,mtags);
    m2id = strcmp(haxes(j).Item2,mtags);
    if ~isequal(old_models(m1id),new_models(m1id)) || ...
        ~isequal(old_models(m2id),new_models(m2id))
        tf(j) = true;
    end
end


% ------------------------------------------------------------------------
function analysisUpdater(~,event) 
%ANALYSISUPDATER Listener callback for getting output from Analysis3DGui
%
% This function is a listener callback which gets invoked during the window
% closing process of ANALYSIS3DGUI.  The listener needs to be added by any
% function which invokes ANALYSIS3DGUI and this function gets the output of
% that GUI and saves/updates all the info in the handles of Registration (on
% successful close)

% Get handles of main GUI & output from closing GUI
[handles,data] = getGuiOutput(event);

% Check output:
if numel(data)==1 && ishandle(data)
    % Figure closed without "OK"
    % Do nothing
    return
end

[hax,loa,marm] = data{:};

% As far as outputs go here, if they have their state variables calculated,
% we can skip them, if not, then we need to recalculate them:
hax_update = cellfun(@isempty,{hax.Axis});
loa_update = cellfun(@isempty,{loa.Point});
marm_update= cellfun(@isempty,{marm.Point1});

% If we have work to do, give some feedback:
if any( [hax_update loa_update marm_update] )
    hl = FigLocker.Lock(handles.figure1);
    hl.setprogress(inf)
end

if any(hax_update)
    hl.settext('Updating helical axes')
    hax(hax_update) = calcHelicalAxes(hax(hax_update),handles.Models);
end

if any(loa_update)
    hl.settext('Updating lines of action')
    loa(loa_update) = calcLinesOfAction(loa(loa_update),handles.Models);
end

if any(marm_update)
    hl.settext('Updating moment arms')
    marm(marm_update) = calcMomentArms(marm(marm_update),handles.Models);
end

% Unlock figure:
if exist('hl','var') == 1
    hl.unlock
end


% Update output variables:
handles.HelicalAxis = hax;
handles.LineOfAction = loa;
handles.MomentArm = marm;

% Update guidata
guidata(handles.figure1,handles)

phaseDisplay(handles);

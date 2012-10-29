function varargout = Solver(varargin)
% SOLVER MATLAB code for Solver.fig
%      SOLVER, by itself, creates a new SOLVER or raises the existing
%      singleton*.
%
%      H = SOLVER returns the handle to a new SOLVER or the handle to
%      the existing singleton*.
%
%      SOLVER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOLVER.M with the given input arguments.
%
%      SOLVER('Property','Value',...) creates a new SOLVER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Solver_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Solver_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Solver

% Last Modified by GUIDE v2.5 18-Oct-2012 14:46:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Solver_OpeningFcn, ...
                   'gui_OutputFcn',  @Solver_OutputFcn, ...
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
end %Solver()


% --- Executes just before Solver is made visible.
function Solver_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Solver (see VARARGIN)

% Prevent re-launch:
if isfield(handles,'initialized')
    return
else
    handles.initialized = 1;
end

% Choose default command line output for Solver
handles.output = hObject;

addpath('Agnes')

% Initialise the tools
initTools(handles,'solver');

% Initialise the solver list:
setgetSolver(handles.SolverPopup);

% Contextual menu:
handles.hscmenu = createContextMenu('cloud');

% Handle the inputs:
if numel(varargin) == 2
    
    % Cloud structures & current phase:
    [models,p] = varargin{:};
    
    % Check form:
    assert(isstruct(models))
        
    % List of names
    tags = {models.Tag}';       
    
    % Number of phases:
    np  = max( cellfun(@numel,{models(:).LoRes}) );
    
    set(handles.ObjectPopup,'String', tags)     % Populate object popup
    set(handles.PhasePopup,'String',(1:np)')    % Populate phase popup
    set(handles.PhasePopup,'Value',p)           % Set phase
    
    % Get the working variables for pose
    %   It's easier to work with the pose information separated from the
    %   model structures while we're in this program, then push back into
    %   models on exit
    [q,x,isdone] = initPose(models);
    
    handles.Models = models;
    handles.q = q;
    handles.x = x;
    handles.isdone = isdone;
    
    % Graphics:
    view(handles.axes1,3)
    axis(handles.axes1,'equal')
    grid(handles.axes1,'on')
    
    % Position gui over the top of calling gui
    positionOver(handles.figure1,gcbf);

    % Now we need to call a graphics update function
    updateDisplay(handles)

else
    warning([mfilename ' should be called with 3 inputs or it will not be functional']);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Solver wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end %Solver_OpeningFcn()


% --- Outputs from this function are returned to the command line.
function varargout = Solver_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end %Solver_OutputFcn()


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end %figure1_ResizeFcn()


% --- Executes on selection change in ObjectPopup.
function ObjectPopup_Callback(hObject, eventdata, handles)
% hObject    handle to ObjectPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ObjectPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ObjectPopup
end %ObjectPopup_Callback()


% --- Executes during object creation, after setting all properties.
function ObjectPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ObjectPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end %ObjectPopup_CreateFcn()


% --- Executes on selection change in PhasePopup.
function PhasePopup_Callback(hObject, eventdata, handles)
% hObject    handle to PhasePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PhasePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PhasePopup
updateDisplay(handles)
end %PhasePopup_Callback()

% --- Executes during object creation, after setting all properties.
function PhasePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PhasePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end %PhasePopup_CreateFcn()


% --- Executes during object creation, after setting all properties.
function EditCloud_button_Callback(hObject, eventdata, handles)
% ... This is all very complicated

% We need a listener which will pop the button up when any of the other
% uitools in the toolbar turn the databrushing toolbar item off:
persistent b li
if ~isa(b,'SelectiveBrush') || ~ishghandle(b.axis)
    b = [];
    %try delete(li), end %#ok<TRYNC>
end

hdisable = [handles.ObjectPopup, handles.PhasePopup];

oj = get(handles.ObjectPopup,'Value');

btnDown = get(hObject,'Value');
if btnDown
    
    % First need to disable changing the phase or object selection, which
    % will stuff things up:
    set(hdisable,'Enable','off')
    
    ojtag = handles.Models(oj).HiRes.Tag;
    htarget = findobj(handles.axes1,'Tag',ojtag);
    b = SelectiveBrush(htarget);
    b.enable
    b.setColor([1 1 1]*0.75)
    
    % Link the "Edit" button to the tool state:
    li = addlistener(b.brushtool,'State',...    
        'PostSet',@(a,b)set(hObject,'Value', isdown(b.AffectedObject) ) );
    li(2) = addlistener(b.brushtool,'State',...
        'PostSet',@(a,b) EditCloud_button_Callback(hObject,[],guidata(hObject)) );
    
else
    % We must delete the listeners first:
    delete(li)
    
    % Then do remaining operations to turn off databrushing:
    set(hdisable,'Enable','on')
    b.disable
    
    data = b.getData;
    clear b
    xyz = data{1};
    
    % If the user pressed "delete", then brushed data will have NaNs in
    % (some?) columns.  Drop these rows if NaNs exist
    nan_rows = any(isnan(xyz),2);
    xyz = xyz(~nan_rows,:);
    
    % Now if data has changed, notify user and push to handles:
    if ~isequal( size(xyz), size(handles.Models(oj).HiRes.xyz) )
        choice = questdlg('These changes will permenantly affect the loaded dataset. Are you sure?',...
            'Modify dataset','Yes','No','No');
        if isequal(choice,'Yes')
            % Modify the model - we use downsample which preserves trace
            % partitioning if it exists
            handles.Models(oj).HiRes = handles.Models(oj).HiRes.downsample(~nan_rows);
            % Graphics & data updates:
            updateDisplay(handles)
            guidata(hObject,handles)
        end
    end
    
end
end %EditCloud_button_Callback()


% --- Executes when selected object is changed in TransformPanel.
function TransformPanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in TransformPanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
end %TransformPanel_SelectionChangeFcn()



% --- Executes on button press in CentreButton.
function CentreButton_Callback(hObject, eventdata, handles)
% hObject    handle to CentreButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get indices
oj = get(handles.ObjectPopup,'Value');
pj = getphase(handles);

% Both a static and a dynamic cloud must exist in order to do matching:
if ~haveClouds(handles.Models(oj),'both')
    return
end

% Cloud shorthand:
dc = handles.Models(oj).LoRes(pj);
s0 = handles.Models(oj).HiRes;

% Current pose:
s1 = s0.transform(handles.q(pj,oj), handles.x(pj,oj,:));

% Create cloud centred on dynamic origin preserving rotation:
deltax = dc.Origin - s1.Origin;
s2 = s1.translate(deltax);

% Record new position:
[q,x] = gettransform(s0,s2,'quat');
handles.q(pj,oj)   = q;
handles.x(pj,oj,:) = x;

% Update:
doUpdates(hObject,handles)
end %CentreButton_Callback()



% --- Executes on button press in AlignCentreButton.
function AlignCentreButton_Callback(hObject, eventdata, handles)
% hObject    handle to AlignCentreButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

oj = get(handles.ObjectPopup,'Value');
pj = getphase(handles);

% Both a static and a dynamic cloud must exist in order to do matching:
if ~haveClouds(handles.Models(oj),'both')
    return
end

% Cloud shorthand:
dc = handles.Models(oj).LoRes(pj);
s0 = handles.Models(oj).HiRes;

% Get orientation:
[q,x] = gettransform(s0,dc,'quat');
handles.q(pj,oj) = q;
handles.x(pj,oj,:) = x;

% Update:
doUpdates(hObject,handles)
end %AlignCentreButton_Callback()



% --- Executes on button press in RevertButton.
function RevertButton_Callback(hObject, eventdata, handles)
% hObject    handle to RevertButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reset cloud:
pj = getphase(handles);
oj = get(handles.ObjectPopup,'Value');

% Static model required for this action:
if ~haveClouds(handles.Models(oj),'static')
    return
end

% Reset orientation:
handles.q(pj,oj) = quaternion([1 0 0 0]);
handles.x(pj,oj,:) = [0 0 0];

% Update:
doUpdates(hObject,handles)
end %RevertButton_Callback()



function AngleEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to AngleEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AngleEditbox as text
%        str2double(get(hObject,'String')) returns contents of AngleEditbox as a double
end %AngleEditbox_Callback()



% --- Executes during object creation, after setting all properties.
function AngleEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AngleEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end %AngleEditbox_CreateFcn()



% --- Executes during object creation, after setting all properties.
function DistanceEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DistanceEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end %DistanceEditbox_CreateFcn()


function DistanceEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to DistanceEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DistanceEditbox as text
%        str2double(get(hObject,'String')) returns contents of DistanceEditbox as a double
end %DistanceEditbox_Callback()



% --- Executes on button press in RotateButton.
function RotateButton_Callback(hObject, eventdata, handles)
% hObject    handle to RotateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get indices
oj = get(handles.ObjectPopup,'Value');
pj = getphase(handles);

% Check we have required static cloud:
if ~haveClouds(handles.Models(oj),'static')
    return
end

% Get the requested rotation:
rang = str2double(get(handles.AngleEditbox,'String'))*pi/180;
raxis = currentComponent(handles,'double');

% Get the current rotation:
sc0 = handles.Models(oj).HiRes;

% Apply present transformations:
sc1 = sc0.transform(handles.q(pj,oj),handles.x(pj,oj,:));

% Then apply the new rotation:
sc2 = sc1.rotate(rang,raxis,sc1.Origin);

% Now find the current q & x and update:
[q,x] = gettransform(sc0,sc2,'quat');
handles.q(pj,oj)   = q; 
handles.x(pj,oj,:) = x;

% Update:
doUpdates(hObject,handles)
end %RotateButton_Callback()


% --- Executes on button press in TranslateButton.
function TranslateButton_Callback(hObject, eventdata, handles)
% hObject    handle to TranslateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get indices:
oj = get(handles.ObjectPopup,'Value');
pj = getphase(handles);

% Check we have required static cloud:
if ~haveClouds(handles.Models(oj),'static')
    return
end

% Retrieve the requested translation:
d = str2double(get(handles.DistanceEditbox,'String'));
v = currentComponent(handles,'double');
deltax = d*v;

% Apply translation:
x0 = handles.x(pj,oj,:);
handles.x(pj,oj,:) = x0(:) + deltax(:);

% Update:
doUpdates(hObject,handles)
end %TranslateButton_Callback()



% --- Executes on selection change in SolverPopup.
function SolverPopup_Callback(hObject, eventdata, handles)
% hObject    handle to SolverPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SolverPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SolverPopup
end %SolverPopup_Callback()


% --- Executes during object creation, after setting all properties.
function SolverPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SolverPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
%end
end %SolverPopup_CreateFcn()



% --- Executes on button press in RefineButton.
function RefineButton_Callback(hObject, eventdata, handles)
% hObject    handle to RefineButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get indices:
oj = get(handles.ObjectPopup,'Value');
pj = getphase(handles);

% Both a static and a dynamic cloud must exist in order to do matching:
if ~haveClouds(handles.Models(oj),'both')
    return
end

% Get clouds:
s0 = handles.Models(oj).HiRes;      % Static cloud
dc = handles.Models(oj).LoRes(pj);  % Dynamic cloud

solver = setgetSolver(handles.SolverPopup);
locker = FigLocker.Lock(handles.figure1);
locker.addprogressbar;
locker.settext(sprintf('Registering %s...',handles.Models(oj).Tag));

% Run the registration which maps s0->dynamic
try
    [q,x] = registerClouds(solver,s0,dc,handles.q(pj,oj),handles.x(pj,oj,:));
catch ME
    locker.unlock;
    rethrow(ME)
end
locker.unlock;

handles.q(pj,oj)   = q;
handles.x(pj,oj,:) = x;

% Update:
doUpdates(hObject,handles)
end %RefineButton_Callback()


% --- Executes on button press in ClearAllButton.
function ClearAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

btn = warncanceldlg('This will erase all calculated registrations. Are you sure you want to continue?',...
    'Erase regisistration');
switch lower( btn )
    case 'ok'
        % Clear all motion
        [handles.Models.q] = deal([]);
        [handles.Models.x] = deal([]);
        
        % Update motion working copies from Models:
        [handles.q,handles.x,handles.isdone] = initPose(handles.Models);
        
        % Update
        doUpdates(hObject,handles)
    case 'cancel'
        return
end
end %ClearAllButton_Callback()


% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure1_CloseRequestFcn(hObject, [], handles)
end %CancelButton_Callback()


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

solver = setgetSolver(handles.SolverPopup);

[q,x] = nanPose(handles.isdone,handles.q,handles.x);

locker = FigLocker.Lock(handles.figure1);
locker.addbuttons('Stop',@regStopFun,'Cancel',@regCancelFun)

[q2,x2,flag] = registerSequence(handles.Models,q,x,solver,@registrationStatus);

locker.unlock;

% Check output - if cancelled, bail out:
if ~flag
    return
end

% Push results back into handles:
handles.isdone = ~q2.isnan;
handles.q = q2;
handles.x = x2;

% Update:
updateDisplay(handles)
guidata(hObject,handles)

% See registerSequence() for configuration of these:
    %-----------------------------------------------
    function regStopFun(~,~)
        global ICPSTOP
        ICPSTOP = 1;
        disp('stop')
    end
    %-----------------------------------------------
    function regCancelFun(~,~)
        global ICPSTOP
        ICPSTOP = -1;
        disp('cancel')
    end
    %-----------------------------------------------
    function registrationStatus(x,s)
        % Update status of registration using locker
        locker.setprogress(x*100)
        locker.settext(s)
    end %registrationStatus()

end %RunButton_Callback()


% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure1_CloseRequestFcn(hObject, [], handles)
end %DoneButton_Callback()


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch hObject
    
    case handles.figure1
        % In this case, user closed using the window X button
        % OUTPUT: []
        handles.output = [];    
        
    case handles.CancelButton
        % User pressed cancel.
        % OUTPUT: []
        handles.output = [];    
        
    case {handles.RunButton, handles.DoneButton}
        % User pressed "Run" or "Done"
        % OUTPUT: Pose information
        
        % For all the positions that have not been calculated numerically
        % (ie, ~isdone), we set their pose info to NaNs:
        [q,x] = nanPose(handles.isdone,handles.q,handles.x);
        
        % If pose for all posisitions are NaN, then there's nothing to
        % keep, so we return empty:
        if all( q.isnan )
            [handles.Models.q] = deal([]);
            [handles.Models.x] = deal([]);
        else
            % If we have any valid data at all, push (q,x) back into the
            % Models structure: 
            for mj = 1:numel(handles.Models)
                handles.Models(mj).q = q(:,mj);
                handles.Models(mj).x = squeeze(x(:,mj,:));
            end
        end
        
        % Return the updated Models structure:
        handles.output = handles.Models;
        
    otherwise
        fprintf(2,'Unknown caller...\n');
        keyboard
end

% Update handles so listener gets up-to-date info:
guidata(hObject,handles)

% First, hide the figure so the user can't mess anything up:
set(handles.figure1,'Visible','off')

% Next, trigger the listener from invoking program, MRIMagic:
set(handles.figure1,'HitTest','off')

% Then when the listener completes, it will return here and delete the
% figure:
delete(handles.figure1)
end %figure1_CloseRequestFcn()






%=========================================================================
%=========================================================================
%=========================================================================


% ------------------------------------------------------------------------
function doUpdates(hObject,handles)
% Run major update functions

% Get/set status flag:
tf = false;
if hObject == handles.RefineButton
    tf = true;
end
oj = get(handles.ObjectPopup,'Value');
pj = getphase(handles);
handles.isdone(pj,oj) = tf;

% Graphics & data updates:
updateDisplay(handles)
guidata(hObject,handles)
end %doUpdates()

% ------------------------------------------------------------------------
function handles = initTools_solver(handles)
% Configure tools and toolbar(s)

hf = handles.figure1;

% Then kill the docking arrow, as we don't want to see it:
set(hf,'DockControls','off')

% Add the standard matlab figure toolbar:
set(hf,'Toolbar','figure')

% Now modify it.  Yair Altman shows on
% http://undocumentedmatlab.com/blog/modifying-default-toolbar-menubar-actions/
% how to find the handles to the tools.  The following two lines show us
% what they all are: 
%
%   hToolbar = findall(gcf,'tag','FigureToolBar');
%   get(findall(hToolbar),'tag')
%
% Which gives:
%     'FigureToolBar'
%     'Plottools.PlottoolsOn'
%     'Plottools.PlottoolsOff'
%     'Annotation.InsertLegend'
%     'Annotation.InsertColorbar'
%     'DataManager.Linking'
%     'Exploration.Brushing'
%     'Exploration.DataCursor'
%     'Exploration.Rotate'
%     'Exploration.Pan'
%     'Exploration.ZoomOut'
%     'Exploration.ZoomIn'
%     'Standard.EditPlot'
%     'Standard.PrintFigure'
%     'Standard.SaveFigure'
%     'Standard.FileOpen'
%     'Standard.NewFigure'

% Now let's find all the ones that we don't want on this GUI:
h = [];
h(end+1) = findall(hf,'tag','Plottools.PlottoolsOn');
h(end+1) = findall(hf,'tag','Plottools.PlottoolsOff');
h(end+1) = findall(hf,'tag','Annotation.InsertColorbar');
h(end+1) = findall(hf,'tag','DataManager.Linking');
h(end+1) = findall(hf,'tag','Standard.EditPlot');
h(end+1) = findall(hf,'tag','Standard.PrintFigure');
h(end+1) = findall(hf,'tag','Standard.SaveFigure');
h(end+1) = findall(hf,'tag','Standard.FileOpen');
h(end+1) = findall(hf,'tag','Standard.NewFigure');

% Destroy all those ones that we don't want:
delete(h);

% Now use our SelectiveBrush class to simply hide the databrushing tool:
SelectiveBrush.HideTool(handles.figure1);

% Now set icon for reset button:
load icons;
set(handles.ClearAllButton,'CData',icons.reset)
end %initTools_solver()


% ------------------------------------------------------------------------
function varargout = setgetSolver(hObj)

if nargout == 0
    % Set up the solver list - in order of preference
    %  These get picked up by register.m
    
    contents = {};
    % LIBICP:
    if 3==exist('icpMex','file') && 3==exist('sparsifyMex','file')
        contents{end+1,1} = 'libicp';
    end
    
    % Jian & Vemuri's Mixtures of Gaussian:
    GMM_matlab = ( exist('gmmreg_L2','file') && exist('gmmreg_L2_costfunc','file') ); % -Matlab implementation
    GMM_bin    = false; % -Compiled binary
    if GMM_matlab || GMM_bin
        contents{end+1,1} = 'Gaussian Mixtures';
    end
    
    % Kjet & Wilm ICP:
    if 2==exist('icp','file')
        contents{end+1,1} = 'Kjer & Wilm ICP';
    end
    
    % UCLA ICP:
    if 2==exist('aditeration','file')
        contents{end+1,1} = 'aditeration';
    end
    
    % Kroon's finite ICP:
    if 2==exist('ICP_finite','file')
        contents{end+1,1} = 'Kroon Finite ICP';
    end    

    set(hObj,'String',contents)
    
elseif nargout == 1
    % Return the currently selected solver
    contents = get(hObj,'String');
    varargout{1} = contents{get(hObj,'Value')};

end
end %setgetSolver()

    
% ------------------------------------------------------------------------
function updateDisplay(handles)


% Shorthand:
hf = handles.figure1;
ha = handles.axes1;

% Retrieve some values:
models = handles.Models;    % Models structure
nm = numel(models);         % Number of models
p = getphase(handles);      % Current phase
clrs = get(ha,'ColorOrder');

% Version specific plot opts:
%   7.12.0 --> 2011a
%   7.13.0 --> 2011b
cprops = {'Marker','.'};            % Default cloud plotting properties
if ~verLessThan('matlab','7.12.0')  % 2011a and later
    % 2011a and later
    cprops = {'Marker','.','MarkerSize',8};
end

% ---------- Delete all objects
delete(get(ha,'Children'));

% Hold on:
hold(ha,'on')

% Draw clouds:
%   - Static or dynamic clouds may exist without the other
for mj = 1:nm
    
    % Static cloud - if one exists and is visible, display it:
    if ~isempty(models(mj).HiRes) && models(mj).HiRes.Visible
        cld = models(mj).HiRes.transform( handles.q(p,mj), handles.x(p,mj,:) );
        hc = cld.plot(ha, 'Color',clrs(mj,:), 'LineStyle','none', cprops{:});
        set(hc,'uicontextmenu',handles.hscmenu)
        cld.plotcs;
    end
    
    % Dynamic cloud - if one exists and is visible, display it:
    if ~isempty(models(mj).LoRes) && models(mj).LoRes(p).Visible
        cld = models(mj).LoRes(p);
        cld.plot(ha,'Color',[1 1 1]*0.1,'Linestyle','none',cprops{:});
        cld.plotcs;
    end
end
axis(ha,'tight')
end %updateDisplay()


% ------------------------------------------------------------------------
function [q,x,isdone] = initPose(models)

% Number of models & phases:
nm = numel(models);
np = max( cellfun(@numel,{models(:).LoRes}) );

static = {models.HiRes};

% Initialise pose matrices
for mj = nm:-1:1
    if isempty(static{mj}) 
        q(1,mj)   = quaternion.NaN;
        x(1,mj,:) = NaN(1,1,3);
    else
        q(1,mj)   = quaternion([1,0,0,0]);  
        x(1,mj,:) = zeros(1,1,3);
    end
end
q = repmat(q(1,:),[np,1]);
x = repmat(x(1,:,:),[np,1,1]);

% Now create blank ISDONE matrix:
isdone = false(np,nm);     % Whether computation has been done or not 

% Now we have our starting matrices which are filled with the pose of the
% static models in their original positions.
%
% From here, we fill those matrices with pose information that has been
% provided as an input, only where it is defined (ie, not NaN).  
%
% First we'll check to make sure that the data provided is ok:

[qin,xin] = models2qx(models,quaternion.NaN,NaN);

if isempty(qin) || isempty(xin)
    return
end
if ~isequal(size(q),size(qin)) || ~isequal(size(x),size(xin))
    error('Input (q,x) are incorrectly sized')
end

% Now assign values:
q( ~isnan(qin) ) = qin( ~isnan(qin) );
x( ~isnan(xin) ) = xin( ~isnan(xin) );
isdone = ~isnan(qin);
end %initPose()


% ------------------------------------------------------------------------
function [qn,xn] = nanPose(isdone,qn,xn)
% NANPOSE 
% For all the positions that have not been calculated numerically
% (ie, ~isdone), we set their pose info to NaNs:
qn(~isdone) = quaternion([NaN NaN NaN NaN]);
xmask = ones(size(isdone));
xmask(~isdone) = NaN;
xn = xn.*repmat(xmask,[1,1,3]);
end %nanPose()


% ------------------------------------------------------------------------
function p = getphase(handles)
hObj = handles.PhasePopup;
contents = cellstr(get(hObj,'String'));
p = str2double(contents{get(hObj,'Value')});
end %getphase()


% ------------------------------------------------------------------------
function tf = haveClouds(model,type)
%HAVECLOUDS Do we have both static & dynamic clouds for this object / phase
% combination so that we can run a match?
tf = true;
msgstr = [];
bMsg = sprintf('No Static or Dynamic clouds loaded for %s',model.Tag);
sMsg = sprintf('No Static cloud loaded for %s',model.Tag);
dMsg = sprintf('No Dynamic clouds loaded for %s',model.Tag);
switch lower( type )
    case 'both'
        if isempty(model.HiRes) && isempty(model.LoRes)
            msgstr = bMsg;
        elseif isempty(model.HiRes)
            msgstr = sMsg;
        elseif isempty(model.LoRes)
            msgstr = dMsg;
        end
    case 'static'
        if isempty(model.HiRes)
            msgstr = sMsg;
        end
    case 'dynamic'
        if isempty(model.LoRes)
            msgstr = dMsg;
        end
    otherwise
        error('Incorrect request');
        
end
% Throw the information popup if required:
if ~isempty(msgstr)
    tf = false;                             % Don't have the requested clouds
    hq = helpdlg(msgstr,'Missing files');   % Display the notification
    uiwait(hq,1.2)                          % Wait for time/destroy
    try delete(hq), end %#ok<TRYNC>         % destroy
    return
end
end %haveClouds()


% ------------------------------------------------------------------------
function cc = currentComponent(handles,opt)
% Get the currently selected ordinate: { x | y | z }
%
% Return in either a string, eg, cc = 'y', or a vector, eg, cc = [0,1,0];
opt = lower(opt);
cc = [get(handles.xRadio,'Value'), get(handles.yRadio,'Value'), get(handles.zRadio,'Value')];

switch opt
    case {'vector','double'}
        return
    case 'string'
        str = [];
        if cc(1)
            str = 'x';
        elseif cc(2)
            str = 'y';
        elseif cc(3)
            str = 'z';
        end
        assert(~isempty(str),'Could not determine the current comoponent');
        cc = str;

end
end %currentComponent()

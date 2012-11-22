function varargout = StateSmoother(varargin)
% STATESMOOTHER MATLAB code for StateSmoother.fig
%      STATESMOOTHER, by itself, creates a new STATESMOOTHER or raises the existing
%      singleton*.
%
%      H = STATESMOOTHER returns the handle to a new STATESMOOTHER or the handle to
%      the existing singleton*.
%
%      STATESMOOTHER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STATESMOOTHER.M with the given input arguments.
%
%      STATESMOOTHER('Property','Value',...) creates a new STATESMOOTHER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StateSmoother_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StateSmoother_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StateSmoother

% Last Modified by GUIDE v2.5 22-Nov-2012 14:15:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StateSmoother_OpeningFcn, ...
                   'gui_OutputFcn',  @StateSmoother_OutputFcn, ...
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


% --- Executes just before StateSmoother is made visible.
function StateSmoother_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StateSmoother (see VARARGIN)

% Choose default command line output for StateSmoother
handles.output = hObject;

% Configure axes:
title(handles.axes1,{'States Vs Phase',' '})
ylabel(handles.axes1,'Angle (deg)')
ylabel(handles.axes2,'Axis Components')
ylabel(handles.axes3,'Position [mm] ')
xlabel(handles.axes3,'Phase')


% GUI *must* be called by "Registration" gui, since its handles contain the
% fields which we will interact with.  If not called by this program,
% StateSmoother will either launch and be non-functional, or will error
% when initialising:
if isempty(gcbf)
    warning([mfilename ' will not be functional when invoked in this way.']) %#ok<WNTAG>
    guidata(hObject,handles)
    return
end

% Store the caller figure handle:
handles.caller = gcbf;
caller_handles = guidata(handles.caller);

% Check we have been invoked by the correct program:
if ~isfield(caller_handles,'Models')
    error(['%s must be called by a subfunction of "Registation.m",',...
        ' since it depends on its handles.'],mfilename)
end

% Get data:
handles.Orig_Models = caller_handles.Models;            % Store initial state of Models
handles.Models = handles.Orig_Models;                   % Then store a copy to work with

handles.Orig_HelicalAxis = caller_handles.HelicalAxis;  % Store inital state
handles.HelicalAxis = caller_handles.HelicalAxis;       % Then store a copy to work with

% Populate listbox
set(handles.listbox1,'String',{handles.Models.Tag})

% Configure the slider:
mn = 0;
mx = max(cellfun(@length,{handles.Models.q}));  % They should all be the same
w  = round(mx/2);   
ss = [0.01 0.1]/( (mx-mn)/10);
set(handles.slider1,'Min',mn)
set(handles.slider1,'Max',mx)
set(handles.slider1,'SliderStep',ss)
set(handles.text_min,'String',num2str(mn))
set(handles.text_max,'String',num2str(mx))

% Set current value - use callbacks to ensure it's within bounds:
set(handles.edit1,'String',num2str(w))
feval(get(handles.edit1,'Callback'),handles.edit1,[])

% Add listeners for refreshing the display:
addlistener(handles.slider1,'Value','PostSet',  @doCalcs);
addlistener(handles.checkbox1,'Value','PostSet',@doCalcs);
addlistener(handles.listbox1,'Value','PostSet', @refreshPlots);


% Update handles structure:
guidata(hObject, handles);

% Run calcs if necessary:
if true % <- condition?
    doCalcs(handles)
else
    % Otherwise just refresh plots:
    refreshPlots(handles)    
end
    
% UIWAIT makes StateSmoother wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = StateSmoother_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.edit1,'String',num2str(get(hObject,'Value')))


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

if get(hObject,'Value')
    state = 'on';
else
    state = 'off';
end

set([handles.slider1,...
    handles.edit1,...
    handles.text_min,...
    handles.text_max],...
    'Enable',state)


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

%set_filter_width(hObject,handles.slider1)

% Check value
v = str2double(get(hObject,'String'));

hs = handles.slider1;

if ~( isnumeric(v) && isfinite(v) )
    v = get(hs,'Value');
    v_old = v;
    
else
    v_old = v;
    mn = get(hs,'Min');
    mx = get(hs,'Max');
    
    % Keep within bounds:
    if v < mn
        v = mn;
    elseif v > mx
        v = mx;
    end
    
end

set(hObject,'String',num2str(v));
set(hs,'Value',v)

% Re-calculate if there has been any value change:
if v ~= v_old
    doCalcs(handles)
end


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
%end


% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure1_CloseRequestFcn(handles.figure1, [], handles)

% --- Executes on button press in OkButton.
function OkButton_Callback(hObject, eventdata, handles)
% hObject    handle to OkButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure1_CloseRequestFcn(handles.figure1, [], handles)



%=========================================================================
%       PRIMARY CALCULATION / DISPLAY FUNCTIONS
%=========================================================================

% ------------------------------------------------------------------------
function doCalcs(varargin)
% DOCALCS is a wrapper for calling functions to re-calculate the
% smoothing, then update the display.
%
% DOCALCS can be called either by a listener, or manually:
% 
%   doCalcs(schema,event)  <- listener call
%   doCalcs(handles)       <- manual call
%

% ----------- Manage inputs -----------
if nargin == 1
    handles = varargin{1};
elseif nargin == 2
    [~,event] = varargin{:};
    handles = guidata(event.AffectedObject);
end


% ----------- Calculate states -----------

SMOOTHING = get(handles.checkbox1,'Value');
if SMOOTHING
    %hl = FigLocker.Lock(handles.figure1);
    %hl.settext('Smoothing data');
    %hl.setprogress(inf);
    %drawnow
    % Smooth the data:
    theta = 1:numel(handles.Orig_Models(1).q);
    span = get(handles.slider1,'Value');
    
    sfun = @(y)smooth(theta,y,span,'rloess');
    handles.Models = handles.Models.smoothpose(sfun);
    
    %hl.unlock;
else
    % Use the raw data:
    handles.Models = handles.Models.clearsmoothing;
end

% ----------- Refresh display -----------

% Update handles:
guidata(handles.figure1,handles)

% Refresh the display:
refreshPlots(handles);

% ----------- Push to main GUI ----------

pushToMain(handles)


% ------------------------------------------------------------------------
function pushToMain(handles)
% PUSHTOMAIN Push the current states to the main program and re-run
% dependent calculations.
%
% PUSHTOMAIN performs the task of updating the main "Registration" GUI with
% the current pose states.  The main implication of this is that
% calculations that depend on these states also need to be re-run.  So the
% process is as follows:
%   1) update helical axis calculations
%   2) push updated Models and HelicalAxis to main gui
%   3) refresh display of main gui

caller_handles = guidata(handles.caller);

% ----------- Re-run dependent calcs -----------

% Helical Axes
hax = caller_handles.HelicalAxis;
hax_update = ~cellfun(@isempty,{hax.Axis});
% If we have work to do, give some feedback:
if any(hax_update)
    %hl.settext('Updating helical axes...')
    fprintf('Updating helical axes...')
    hax(hax_update) = calcHelicalAxes(hax(hax_update),handles.Models);
    disp('done!')
end
caller_handles.HelicalAxis = hax;

% Moment Arms
%  - or are these fast enough to do on demand? (ie, in phaseDisplay)

% ----------- Push data -----------

% As the last operation, push the updated Models back to the invoking
% program:
caller_handles.Models = handles.Models;
guidata(handles.caller,caller_handles);

% ----------- Refresh main display -----------

% Run callback to refresh the view:
obj = caller_handles.PhaseSlider;
cbk = get(obj,'Callback');
cbk(obj,[]);



% ------------------------------------------------------------------------
function refreshPlots(varargin)
% REFRESHPLOTS Refresh the plots with the current pose states, given the
% currently selected model

% Manage inputs
if nargin == 1
    handles = varargin{1};
elseif nargin == 2
    [~,event] = varargin{:};
    handles = guidata(event.AffectedObject);
end

% Clear axes:
axset = [handles.axes1, handles.axes2, handles.axes3];
for axj = axset
    hold(axj,'on')
    cla(axj)
end

% Get the currently selected bone:
b = handles.Models(get(handles.listbox1,'Value'));

% Configure colours:
if b.smoothed
    % If plotting smoothed states, plot raw in grey & smooth in colour
    rclrs = lines(4);%ones(4,3)*0.5;
    rwid = 0.5;
    lsty = '--';
    sclrs = lines(4);
    swid  = 1.5;
    
else
    % If only plotting raw states, plot raw in colour
    rclrs = lines(4);
    rwid  = 1.5;
    lsty = '-';
end

% Get raw states:
R2D = 180/pi;
[r_ang,r_axs] = b.qraw.angleaxis;
r_ang = r_ang*R2D;
r_xyz = b.xraw;

% Plot raw states:
rprops = {'LineStyle',lsty,'LineWidth',rwid};
plot(axset(1),r_ang,     'color',rclrs(1,:),rprops{:})

plot(axset(2),r_axs(1,:),'color',rclrs(1,:),rprops{:})
plot(axset(2),r_axs(2,:),'color',rclrs(2,:),rprops{:})
plot(axset(2),r_axs(3,:),'color',rclrs(3,:),rprops{:})
ylim(axset(2),[-1 1])

plot(axset(3),r_xyz(:,1),'color',rclrs(1,:),rprops{:})
plot(axset(3),r_xyz(:,2),'color',rclrs(2,:),rprops{:})
plot(axset(3),r_xyz(:,3),'color',rclrs(3,:),rprops{:})

% Plot smoothed states if necessary:
if b.smoothed
    % Plot smoothed states over the top
    [s_ang,s_axs] = b.qsmooth.angleaxis;
    s_ang = s_ang*R2D;
    s_xyz = b.xsmooth;
    
    props = {'LineWidth',swid};
    plot(axset(1),s_ang,'color',sclrs(1,:),props{:})
    
    plot(axset(2),s_axs(1,:),'color',sclrs(1,:),props{:})
    plot(axset(2),s_axs(2,:),'color',sclrs(2,:),props{:})
    plot(axset(2),s_axs(3,:),'color',sclrs(3,:),props{:})
    ylim(axset(2),[-1 1])
    
    plot(axset(3),s_xyz(:,1),'color',sclrs(1,:),props{:})
    plot(axset(3),s_xyz(:,2),'color',sclrs(2,:),props{:})
    plot(axset(3),s_xyz(:,3),'color',sclrs(3,:),props{:})
end



%=========================================================================
%       HELPER FUNCTIONS
%=========================================================================


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


switch gcbo
    
    case {handles.figure1, handles.OkButton}
        % figure1  - User closed with window "X"
        % OkButton - User closed with "Ok"
        
        % Leave everythign as it is
    case handles.CancelButton
        % CancelButton - User closed with "Cancel"
        
        % Reset Models & HelicalAxis
        handles.Models = handles.Orig_Models;
        handles.HelicaAxis = handles.Orig_HelicalAxis;
                
        % Push back to main gui:
        pushToMain(handles)
        
    otherwise
        error('Unhandled close method')
        
end
        
        



% Hint: delete(hObject) closes the figure
delete(hObject);

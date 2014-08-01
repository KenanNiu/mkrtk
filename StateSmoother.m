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

% Last Modified by GUIDE v2.5 27-Nov-2012 12:01:00

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
ylabel(handles.axes1,'Euler Angles (deg)')
ylabel(handles.axes2,'Position [mm] ')
xlabel(handles.axes2,'Phase')

% Add the radio boxes:
handles = newFilterRadios(handles);
handles = newAbscissaRadios(handles);

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

% Copy relevant data from main program to handles of StateSmoother:
handles.Orig.Models = caller_handles.Models;            % Store initial state of Models
handles.Models      = handles.Orig.Models;              % Then store a copy to work with

handles.Orig.HelicalAxis = caller_handles.HelicalAxis;  % Store inital state of Axes
handles.HelicalAxis      = handles.Orig.HelicalAxis;    % Then store a copy to work with

if isfield(caller_handles,'pose_filter') && ~isempty(caller_handles.pose_filter)
    handles.Orig.pose_filter = caller_handles.pose_filter; % Store initial state
    handles.pose_filter      = handles.Orig.pose_filter;   % Then store a copy to work with
else
    handles.Orig.pose_filter.enabled = false;
    handles.pose_filter.enabled      = false;
end

% Populate listbox
set(handles.Listbox_Models,'String',{handles.Models.Tag})

% Configure the slider:
mn = 0;
mx = max(cellfun(@length,{handles.Models.q}));  % They should all be the same
  
ss = [0.01 0.1]/( (mx-mn)/10);
set(handles.Slider_Width,'Min',mn)
set(handles.Slider_Width,'Max',mx)
set(handles.Slider_Width,'SliderStep',ss)
set(handles.Text_Min,'String',num2str(mn))
set(handles.Text_Max,'String',num2str(mx))

positionOver(handles.figure1,caller_handles.figure1)

% Update handles structure:
guidata(hObject, handles);


% Apply filter settings from pose_filter:
applyFilterSettings(handles)

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
function Slider_Width_Callback(hObject, eventdata, handles)
% hObject    handle to Slider_Width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Update edit box:
set(handles.Edit_Width,'String',num2str(get(hObject,'Value')))

% Re-calculate
doCalcs(handles)


% --- Executes during object creation, after setting all properties.
function Slider_Width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slider_Width (see GCBO)
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


% ---- Enable/disable uicontrols ---- %
if get(hObject,'Value')
    state = 'on';
else
    state = 'off';
end

% Get handles to radio buttons:
hr = findall(handles.figure1,'Type','uicontrol','style','radiobutton','-regexp','Tag','radio_');

% Now disable/enable all the dependent uicontrols:
set([handles.Text_FilterWidth
    handles.Slider_Width
    handles.Edit_Width
    handles.Text_Min
    handles.Text_Max
    handles.Text_Filter
    hr],...
    'Enable',state)

% ---- Update / refresh ---- %
doCalcs(handles)


function Edit_Width_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_Width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edit_Width as text
%        str2double(get(hObject,'String')) returns contents of Edit_Width as a double

%set_filter_width(hObject,handles.Slider_Width)

% Check value
v = str2double(get(hObject,'String'));

hs = handles.Slider_Width;

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
function Edit_Width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_Width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Listbox_Models.
function Listbox_Models_Callback(hObject, eventdata, handles)
% hObject    handle to Listbox_Models (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Listbox_Models contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Listbox_Models

doCalcs(handles)

% --- Executes during object creation, after setting all properties.
function Listbox_Models_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Listbox_Models (see GCBO)
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


function AbscissaRadios_Callback(hObject, eventdata, handles)
% Always keep one radio button current:
set(hObject,'Value',1)
hp = get(hObject,'Parent');
hr = findall(hp,'Type','uicontrol','Style','radio','-regexp','tag','radio_abscissa_');
set(hr(hr~=hObject),'Value',0)
drawnow

% Run calcs
doCalcs(handles)


function FilterRadios_Callback(hObject, eventdata, handles)
% Always keep one radio button current:
set(hObject,'Value',1)
hp = get(hObject,'Parent');
hr = findall(hp,'Type','uicontrol','Style','radio','-regexp','tag','radio_filter_');
set(hr(hr~=hObject),'Value',0)
drawnow

% Run calcs
doCalcs(handles)


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
        
        % Reset Models & HelicalAxis & pose_filter
        handles.Models = handles.Orig.Models;
        handles.HelicaAxis = handles.Orig.HelicalAxis;
        handles.pose_filter = handles.Orig.pose_filter;
                
        % Push back to main gui:
        pushToMain(handles)
        
    otherwise
        error('Unhandled close method')
        
end

% Hint: delete(hObject) closes the figure
delete(hObject);


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
    % Get filter states from GUI:
    [theta,slider_val,sfun] = guiFilterStates(handles);
    
    % Smooth models
    handles.Models = handles.Models.smoothpose(sfun);
    
    % Store pose_filter states:
    handles.pose_filter.enabled = true;
    handles.pose_filter.theta = theta;
    handles.pose_filter.slider_val  = slider_val;
    handles.pose_filter.fun   = sfun;
    
    %hl.unlock;
else
    % Use the raw data:
    handles.Models = handles.Models.clearsmoothing;
    
    % Store filter props:
    handles.pose_filter.enabled = false;
end

% ----------- Refresh display -----------

% Update handles:
guidata(handles.figure1,handles)

% Refresh the display:
refreshPlots(handles);

% ----------- Push to main GUI ----------

pushToMain(handles)


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
axset = [handles.axes1, handles.axes2];
for axj = axset
    hold(axj,'on')
    delete(findall(axj,'type','line'));
end

% Get the currently selected bone:
b = handles.Models(get(handles.Listbox_Models,'Value'));

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
r_eul = b.qraw.eulerangles('xyz')*R2D;
r_xyz = b.xraw;

% Plot raw states:
rprops = {'LineStyle',lsty,'LineWidth',rwid};

hi = plot(axset(1),r_eul(1,:),'color',rclrs(1,:),rprops{:});
hj = plot(axset(1),r_eul(2,:),'color',rclrs(2,:),rprops{:});
hk = plot(axset(1),r_eul(3,:),'color',rclrs(3,:),rprops{:});

hx = plot(axset(2),r_xyz(:,1),'color',rclrs(1,:),rprops{:});
hy = plot(axset(2),r_xyz(:,2),'color',rclrs(2,:),rprops{:});
hz = plot(axset(2),r_xyz(:,3),'color',rclrs(3,:),rprops{:});

% Plot smoothed states if necessary:
if b.smoothed
    set([hi,hj,hk, hx,hy,hz],'HandleVisibility','off')
    
    % Plot smoothed states over the top
    s_eul = b.qsmooth.eulerangles('123')*R2D;
    s_xyz = b.xsmooth;
    
    props = {'LineWidth',swid};
    
    hi = plot(axset(1),s_eul(1,:),'color',sclrs(1,:),props{:});
    hj = plot(axset(1),s_eul(2,:),'color',sclrs(2,:),props{:});
    hk = plot(axset(1),s_eul(3,:),'color',sclrs(3,:),props{:});
    
    hx = plot(axset(2),s_xyz(:,1),'color',sclrs(1,:),props{:});
    hy = plot(axset(2),s_xyz(:,2),'color',sclrs(2,:),props{:});
    hz = plot(axset(2),s_xyz(:,3),'color',sclrs(3,:),props{:});
end

set(hi,'DisplayName','roll')
set(hj,'DisplayName','pitch')
set(hk,'DisplayName','yaw')

set(hx,'DisplayName','x')
set(hy,'DisplayName','y')
set(hz,'DisplayName','z')

% Configure legends
legend(axset(1),'show');
legend(axset(2),'show');

drawnow

%=========================================================================
%       HELPER FUNCTIONS
%=========================================================================


% ------------------------------------------------------------------------
function handles = addToHandles(handles,hobjs)
% Add the uicontrol handles in HOBJS to HANDLES:
for j = 1:numel(hobjs)
    tag = get(hobjs(j),'Tag');
    handles.(tag) = hobjs(j);
end


% ------------------------------------------------------------------------
function applyFilterSettings(handles)
% Apply filter settings to the gui
%

% First we need to update handles, otherwise we can't run the callbacks we
% need to (since GUI is still building):
guidata(handles.figure1,handles)

fprops = handles.pose_filter;

% Enabled / disabled:
set(handles.checkbox1,'Value',fprops.enabled)

% Abscissa
if isfield(fprops,'theta') && fprops.enabled
    if all(diff(fprops.theta) == 1)
        aObj = handles.radio_abscissa_phaseId;
    else
        aObj = handles.radio_abscissa_angle;
    end
    set(aObj,'Value',1)
else
    set(handles.radio_abscissa_phaseId,'Value',1)
end

% Smoothing function:
regtag = 'radio_filter_';
rObjs = findall(handles.figure1,...
    'Type','uicontrol',...
    'Style','radiobutton',...
    '-regexp','Tag',regtag);
fnames = regexprep(get(rObjs,'Tag'),regtag,'');
if isfield(fprops,'fun') && fprops.enabled
    funstr = char(fprops.fun);
    n = 0;
    i = 1;
    for j = 1:numel(fnames)
        [i1,i2] = regexp(funstr,fnames{j});
        if (i2-i1) > n
            n = i2-i1;
            i = j;
        end
    end
    set(rObjs(i),'Value',1)
else
    set(rObjs(1),'Value',1)
end

% Smoothing value
sld = handles.Slider_Width;
if isfield(fprops,'slider_val') && fprops.enabled
    sval = fprops.slider_val;    
else
    sval = round(mean([get(sld,'Min'), get(sld,'Max')]));
end
set(handles.Edit_Width,'String',num2str(sval))
set(handles.Slider_Width,'Value',sval)

% Then run the callback for the checkbox, which will run the updates:
runCallback(handles.checkbox1)


% ------------------------------------------------------------------------
function handles = newAbscissaRadios(handles)
% Create new radio buttons for 

hHang = handles.Text_Abscissa;
hparent = get(hHang,'Parent');

% Base hang point:
psn = get(hHang,'Position');

% Set width & height:
w = 150;
h = 20;
vsep = 1;
indent = 20;

% Set initial position & increment:
% Set initial position & increment:
psn(1:2) = [psn(1)+indent psn(2)-h-vsep*2];
psn(3:4) = [w,h];
increment = @(p)[p(1) p(2)-vsep-h psn(3) psn(4)];

% Create a convenient callback generator using an anonymous fcn:
%   (Makes a callback in the standard format used by GUIDE)
cbk = eval(['@(hObject,eventdata)' mfilename ...
    '(''AbscissaRadios_Callback'',hObject,eventdata,guidata(hObject))']);
props = {'Style','radiobutton','Callback',cbk};

h = [];

h(end+1) = uicontrol(hparent,props{:},'String','Phase ID','tag','radio_abscissa_phaseId','Position',psn);
psn = increment(psn);
h(end+1) = uicontrol(hparent,props{:},'String','Estimated joint angle','tag','radio_abscissa_angle','Position',psn);

handles = addToHandles(handles,h);


% ------------------------------------------------------------------------
function handles = newFilterRadios(handles)
% Create radio buttons for choosing filter

hHang = handles.Text_Filter;
hparent = get(hHang,'Parent');

% Base hang point:
psn = get(hHang,'Position');

% Set width & height:
w = 150;
h = 20; 
vsep = 1;       % [pix]
indent = 20;    % [pix]

% Set initial position & increment:
psn(1:2) = [psn(1)+indent psn(2)-h-vsep*2];
psn(3:4) = [w,h];
increment = @(p)[p(1) p(2)-vsep-h psn(3) psn(4)];

% Create a convenient callback generator using an anonymous fcn:
%   (Makes a callback in the standard format used by GUIDE)
cbk = eval(['@(hObject,eventdata)' mfilename ...
    '(''FilterRadios_Callback'',hObject,eventdata,guidata(hObject))']);
props = {'Style','radiobutton','Callback',cbk};
h = [];

% Add their names to handles??
h(end+1) = uicontrol(hparent,props{:},'String','Moving average','tag','radio_filter_moving','Position',psn);

if exist('csaps','file') == 2
    psn = increment(psn);
    h(end+1) = uicontrol(hparent,props{:},'String','Cubic smoothing spline','tag','radio_filter_csaps','Position',psn);
end

if exist('smooth','file') == 2
    psn = increment(psn);
    h(end+1) = uicontrol(hparent,props{:},'String','Loess','tag','radio_filter_loess','Position',psn);
    psn = increment(psn);
    h(end+1) = uicontrol(hparent,props{:},'String','Lowess','tag','radio_filter_lowess','Position',psn);
    psn = increment(psn);
    h(end+1) = uicontrol(hparent,props{:},'String','Robust loess','tag','radio_filter_rloess','Position',psn);
    psn = increment(psn);
    h(end+1) = uicontrol(hparent,props{:},'String','Robust lowess','tag','radio_filter_rlowess','Position',psn);
end

handles = addToHandles(handles,h);


% ------------------------------------------------------------------------
function [theta,slider_val,sfun] = guiFilterStates(handles)
% Get the filter states from the GUI

% Get filter span:
slider_val = get(handles.Slider_Width,'Value');

% Get the independent variable - theta:
theta = [];
if get(handles.radio_abscissa_angle,'Value')    % If 'angle' is selected
    if isfield(handles,'HelicalAxis') && ...       % If HelicalAxis field exists
        ~isempty(handles.HelicalAxis)           %   ...and is populated
    theta = joint_dt_proxy(handles.Models);     % Then get joint angle from that
    else
        set(handles.radio_abscissa_angle,  'Value',0)   % Switch to phase
        set(handles.radio_abscissa_phaseId,'Value',1)   %
    end
end
% But if the above didn't yield anything, either 'phase' is selected, or
% there are not Helical Axes to work from
if isempty(theta) %get(handles.radio_abscissa_phaseId,'Value')
    theta = 1:numel(handles.Orig.Models(1).q);
end


% For filters that use 0->1, they tend to fail on 0 and 1, so limit the
% value to the range [eps 1-eps]:
num2frac = @(v)min([max([v/get(handles.Slider_Width,'Max'), eps]), 1-eps]);

% Radio button tags will have the form: 
%   'radio_filter_lowess'
% so we just grab out the name at the end for the switch statement below.
searchstr = 'radio_filter_';
hr = findall(handles.figure1,'-regexp','Tag',searchstr,'Value',1);
funstr = regexprep(get(hr,'tag'),searchstr,'');

span = slider_val;
switch funstr
    case 'moving'
        % This fails if span < 1, so limit it, but don't bother about
        % feeding it back to the GUI.
        if span < 1
            span = 1;
        end
        sfun = @(y)smooth(theta,y,double(span),'moving');
        %window = ones(span,1)/span;
        %sfun = @(y)convn(y,window,'same');
        
    case 'csaps'
        f = 1-num2frac(span);   % Fraction, in the range: 0 -> 1
        b = f*2-1;              % Bilateral, in the  range: -1 -> 1
        h = mean(diff(theta));  % Average datasite spacacing (see doc csaps)
        e = 3;                  % Exponential factor
        den = 10^(b*e)*6;       % Denominator term (see doc csaps)
        p = 1/(1 + h^3/den);    % Smoothing parameter (see doc csaps)
        span = p;               % Now this is our new equivalent term for csaps
        sfun = @(y)csaps(theta,y,span,theta);
        
    case 'loess'
        span = num2frac(span);
        sfun = @(y)smooth(theta,y,span,'loess');
        
    case 'lowess'
        span = num2frac(span);
        sfun = @(y)smooth(theta,y,span,'lowess');
        
    case 'rloess'
        span = num2frac(span);
        sfun = @(y)smooth(theta,y,span,'rloess');
        
    case 'rlowess'
        span = num2frac(span);
        sfun = @(y)smooth(theta,y,span,'rlowess');
end



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
if ~isempty(hax)
    hax_update = ~cellfun(@isempty,{hax.Axis});
    % If we have work to do, give some feedback:
    if any(hax_update)
        %hl.settext('Updating helical axes...')
        fprintf('Updating helical axes...')
        hax(hax_update) = calcHelicalAxes(hax(hax_update),handles.Models);
        disp('done!')
    end
    caller_handles.HelicalAxis = hax;
end
% Moment Arms
%  - or are these fast enough to do on demand? (ie, in phaseDisplay)

% ----------- Push data -----------

% As the last operation, push the updated Models back to the invoking
% program:
caller_handles.Models = handles.Models;
% Also push pose_filter states:
caller_handles.pose_filter = handles.pose_filter;
% Update:
guidata(handles.caller,caller_handles);

% ----------- Refresh main display -----------

% Run callback to refresh the view:
obj = caller_handles.PhaseSlider;
cbk = get(obj,'Callback');
cbk(obj,[]);

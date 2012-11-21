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

% Last Modified by GUIDE v2.5 22-Nov-2012 09:08:39

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

% Store the caller figure handle:
handles.caller = gcbf;

% Configure axes:
ylabel(handles.axes1,'Angle (deg)')
ylabel(handles.axes2,'Axis Components')
ylabel(handles.axes3,'Position [mm] ')
xlabel(handles.axes3,'Phase')

% Check for call with no inputs - gui not functional, but circumvent error:
if numel(varargin) == 0
    guidata(hObject,handles)
    return
end    

% Get inputs:
handles.Models_orig = varargin{1};      % Store initial state of Models
handles.Models = handles.Models_orig;   % Then store a copy to work with

% Populate listbox
set(handles.listbox1,'String',{handles.Models.Tag})


% Add listeners for refreshing the display:
addlistener(handles.slider1,'Value','PostSet',@doCalcs);
addlistener(handles.listbox1,'Value','PostSet',@doCalcs);
addlistener(handles.checkbox1,'Value','PostSet',@doCalcs);


% Update handles structure:
guidata(hObject, handles);

% Run calcs if necessary:
if true
    disp('fix condition')
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
disp('Cancel - no action defined yet')

% --- Executes on button press in OkButton.
function OkButton_Callback(hObject, eventdata, handles)
% hObject    handle to OkButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Ok - no action defined yet')


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
    disp('**** lock fig? ****')
    % Smooth the data:
    %handles.Models = smooth_models(handles.Models_orig,x,span);
    disp('**** fix this ****')
    handles.Models = handles.Models_orig;
else
    % Use the raw data:
    handles.Models = handles.Models_orig;
end

% ----------- Refresh display -----------

% Update handles:
guidata(handles.figure1,handles)

% Refresh the display:
refreshPlots(handles);


% ----------- Push data -----------


% As the last operation, push the updated Models back to the invoking
% program:
%caller_handles = guidata(handles.caller);
%caller_handles.Models = handles.Models;
%guidata(handles.caller,caller_handles);




% ------------------------------------------------------------------------
function refreshPlots(handles)
% REFRESHPLOTS Refresh the plots with the current pose states, given the
% currently selected model

axset = [handles.axes1, handles.axes2, handles.axes3];
for axj = axset
    hold(axj,'on')
    cla(axj)
end

% Get the currently selected bone:
b = handles.Models(get(handles.listbox1,'Value'));

% Configure colours:
if b.smoothed
    rclrs = zeros(4,3);
    % plot raw states
    % Overlay smoothed states
else
    rclrs = lines(4);
end

% Get raw states:
R2D = 180/pi;
[r_ang,r_axs] = b.qraw.angleaxis;
r_ang = r_ang*R2D;
r_xyz = b.xraw;

% Plot raw states:
plot(axset(1),r_ang,'color',rclrs(1,:))

plot(axset(2),r_axs(1,:),'color',rclrs(1,:))
plot(axset(2),r_axs(2,:),'color',rclrs(2,:))
plot(axset(2),r_axs(3,:),'color',rclrs(3,:))
ylim(axset(2),[-1 1])

plot(axset(3),r_xyz(1,:),'color',rclrs(1,:))
plot(axset(3),r_xyz(2,:),'color',rclrs(2,:))
plot(axset(3),r_xyz(3,:),'color',rclrs(3,:))

% Plot smoothed states if necessary:
if b.smoothed
    % plot smoothed states
    [s_ang,s_axs] = b.qsmooth.angleaxis;
end










%=========================================================================



% ------------------------------------------------------------------------
function mdls = smooth_models(mdls,theta,span)

n = numel(mdls);
for j = 1:n     % ==> upgrate to parfor
    if ~isempty(mdls(j).q) && ~all(isnan(mdls(j).q))
        qj = mdls(j).qraw.unwrap;   %\_ raw data
        %qj = mdls(j).qraw;
        xj = mdls(j).xraw;          %/
        
        sfun = @(y)smooth(theta,y,span,'rloess');
        
        qj = qj.smooth(sfun,1);
        for c = 1:3
            xj(:,c) = sfun(xj(:,c));
        end
        
        mdls(j).q = qj;
        mdls(j).x = xj;
    end
end

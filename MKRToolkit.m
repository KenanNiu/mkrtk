function varargout = MKRToolkit(varargin)
% MKRTOOLKIT MATLAB code for MKRToolkit.fig
%      MKRTOOLKIT, by itself, creates a new MKRTOOLKIT or raises the existing
%      singleton*.
%
%      H = MKRTOOLKIT returns the handle to a new MKRTOOLKIT or the handle to
%      the existing singleton*.
%
%      MKRTOOLKIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MKRTOOLKIT.M with the given input arguments.
%
%      MKRTOOLKIT('Property','Value',...) creates a new MKRTOOLKIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MKRToolkit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MKRToolkit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MKRToolkit

% Last Modified by GUIDE v2.5 24-Oct-2012 14:49:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MKRToolkit_OpeningFcn, ...
                   'gui_OutputFcn',  @MKRToolkit_OutputFcn, ...
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


% --- Executes just before MKRToolkit is made visible.
function MKRToolkit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MKRToolkit (see VARARGIN)

% Choose default command line output for MKRToolkit
handles.output = hObject;

% Position nicely
movegui(handles.figure1,'center')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MKRToolkit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MKRToolkit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SegmentationPushbutton.
function SegmentationPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to SegmentationPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Segmentation

% --- Executes on button press in RegistrationPushbutton.
function RegistrationPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to RegistrationPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Registration

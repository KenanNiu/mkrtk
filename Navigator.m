function varargout = Navigator(varargin)
% NAVIGATOR MATLAB code for Navigator.fig
%      NAVIGATOR, by itself, creates a new NAVIGATOR or raises the existing
%      singleton*.
%
%      H = NAVIGATOR returns the handle to a new NAVIGATOR or the handle to
%      the existing singleton*.
%
%      NAVIGATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NAVIGATOR.M with the given input arguments.
%
%      NAVIGATOR('Property','Value',...) creates a new NAVIGATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Navigator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Navigator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Navigator

% Last Modified by GUIDE v2.5 27-Mar-2012 11:08:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Navigator_OpeningFcn, ...
    'gui_OutputFcn',  @Navigator_OutputFcn, ...
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


% --- Executes just before Navigator is made visible.
function Navigator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Navigator (see VARARGIN)

% Prevent re-launch:
if isfield(handles,'initialized')
    return
else
    handles.initialized = 1;
end

% Choose default command line output for Navigator
handles.output = hObject;
handles.caller = gcbf;

initTools(handles.figure1)

% Get java frame - use timer object to do this once figure is visible
handles.jframe = [];
t = timer('TimerFcn',@setjframe,...     % Set handles.jframe
    'StartDelay',0.1,...                % Ensures figure will be visible
    'ExecutionMode','singleShot',...    % Fire just once
    'UserData',handles.figure1);        % Pass it the figure handle


% Update Display:
view(handles.axes1,3)
handles.gap = 0.5;  % gap between bars
handles = updateDisplay(handles);

% Add listener which will update this display whenever the main program
% updates its display.  We hang it on the userdata of the image in the
% main program since this is the method UPDATESLICE provides for triggering
% listeners. We also store the handle to the listener so that it can be
% deleted when this figure is closed:
him = findobj(handles.caller,'Type','image');
handles.listener = addlistener(him,'UserData','PostSet',@updateDisplay);
handles.listener.CallbackTarget = handles.figure1;  % Store this

% Update handles structure
guidata(hObject, handles);

% Now fire the timer which sets handles.jframe then self-destructs:
start(t)

% UIWAIT makes Navigator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Navigator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Delete the listener:
delete(handles.listener);

% Delete the figure:
delete(hObject);


% ------------------------------------------------------------------------
function initTools(hf)
% Configure toolbars & other tools:


% Firstly add the standard matlab figure toolbar:
set(hf,'Toolbar','figure')

% Then kill the docking arrow, as we don't want to see it:
set(hf,'DockControls','off')

% Now modify the standard toolbar.  Yair Altman shows on
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

set(findall(hf,'tag','Standard.EditPlot'),'Visible','off');

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

% Remove separator:
set(findall(hf,'Tag','Exploration.ZoomIn'),'Separator','off')


% ------------------------------------------------------------------------
function setjframe(varargin)
%SETJFRAME Timer function to set handles.jframe
% The reason we need to use a timer and this timer function rather than a
% listener is because the 'Visible' property of a figure gets set before
% the figure actually becomes visible.  So we simply wait for a specified
% interval, 50-100ms, then do the job.
[tobj,~] = varargin{:};     % { tobj, tdata }
hf = get(tobj,'UserData');

% Get java handle & update:
handles = guidata(hf);
handles.jframe = getjframe(hf);
guidata(hf,handles)

% Self-destruct:
stop(tobj)
delete(tobj);


% ------------------------------------------------------------------------
function refreshControls(handles)

% Delete old controls

% Create new controls

% ------------------------------------------------------------------------
function [nslices,nphases,nobjects,uqclrs,slices,phases,objects,clrs] = ...
    getCurrentStats(hMainGui)
handles = guidata(hMainGui);          % Main GUI handles
traces = handles.traces;              % Retrieve traces
clrs = ({handles.traces.Color});      % Colours list

nslices = size(handles.Images.info,1); % Number of slices
nphases = size(handles.Images.info,2); % Number of phases

slices = [traces.Slice]';           % List of traced slices
phases = [traces.Phase]';           % List of traced phases

% Now we need to convert the colours to proper rgb colours:
clrs = ColourSpec('ToRGB',clrs);

% Now find number of objects & oject numbers from the colours:
uqclrs = unique(clrs,'rows');
nobjects = size(uqclrs,1);                % Number of objects
objects = zeros(size(slices));            % Empty object list
for oid = 1:nobjects
    C = repmat(uqclrs(oid,:),numel(slices),1);  % Temporary colour matrix
    objects(all(clrs==C,2)) = oid;              % Fill with object ids
end


% ------------------------------------------------------------------------
function handles = updateDisplay(varargin)
%UPDATEDISPLAY Function which re-draws all items on axes1
%
% Usage:
% ------
%   updateDisplay(handles)          % For a manual call
%   updateDisplay(schema,event)     % Call from a listener
%   updateDisplay(hObject,event)    % Call from listener with 'CallbackTarget' set
%
%
% If the listener property 'CallbackTarget' is set when creating the
% listener, that value comes in as the first input instead of getting a
% schema structure.
%
% NOTE:
% ------
%   On every call, this function destroys and re-creates all the cubes.
%   This is a slow process and could be sped up by doing the following:
%   - Store info from last time this function was run, namely
%       - Number of slices
%       - Number of phases
%       - Number of objects
%       - Number of ROIs
%   - If none of the above have changed, we don't need to do anything
%   - If number of ROIs has changed, need to colour/uncolour cubes
%   - If any of the others have changed, we do a full re-build
%
%  This will speed things up significantly

if nargin == 1
    % Manual call: updateDisplay(handles)
    handles = varargin{1};
elseif nargin == 2
    % Listener call: updateDisplay(figure,event)
    handles = guidata(varargin{1});
end


% The different colours in the plot identify different bones/objects.

OBJSEP = 5;
zvalue = @(oj)(oj-1)*OBJSEP + 1;

% Some shorthand
hf = handles.figure1;
hax = handles.axes1;

% Get current stats:
[ns,np,no,uqclrs,s,p,o,clrs] = getCurrentStats(handles.caller);

% Add/refresh uicontrols:
refreshControls(handles);

% Clear axes:
delete(findobj(handles.axes1,'Type','patch'));

% IsTraced matrix: slices in rows, phases in columns, objects in layers
IT = false(ns,np,no);
idx = sub2ind(size(IT),s,p,o);
IT(idx) = true;

% Plot data:
zc = zvalue(o);
h1 = drawCube(hax,[s(:),p(:),zc(:)],o(:),clrs,clrs/1.5);

% NotTraced matrix, and associated indices:
NT = ~IT;
[s,p,o] = ind2sub(size(NT),find(NT));
clrs = uqclrs(o,:);

% Plot data
zc = zvalue(o);
h0 = drawCube(hax,[s(:),p(:),zc(:)],o(:),ones(numel(s),3)*.95,clrs/1.5);


% Format axes & labels
xlabel(hax,'Slice')
xtick = 1:ns;
set(hax,'XTick',xtick-0.5)
set(hax,'XTickLabel',xtick)

ylabel(hax,'Phase')
ytick = 1:np;
set(hax,'YTick',ytick-0.5)
set(hax,'YTickLabel',ytick)


zlabel(hax,'Object')
ztick = zvalue(1:no);
set(hax,'ZTick',ztick-0.5)
set(hax,'ZTickLabel',1:no)


axis(hax,'equal','tight')

% Position gui
if ~isequal('on',get(hf,'Visible'))
    movegui(hf,'center')
end

% Add callbacks:
set([h1; h0],'ButtonDownFcn',@cubeBDF)
%set(hf,'WindowButtonDownFcn',@selectWBDF)


% ------------------------------------------------------------------------
function hcube = drawCube(hax,ctrxyz,objid,faceClr,edgeClr)
vertices = [...
    0 0 0
    1 0 0
    1 1 0
    0 1 0
    0 0 1
    1 0 1
    1 1 1
    0 1 1 ];
vertices = vertices.*.8+0.1;
faces = [...
    1 2 6 5
    2 3 7 6
    3 4 8 7
    4 1 5 8
    1 2 3 4
    5 6 7 8];

nc = size(ctrxyz,1);
hcube = zeros(nc,1);
for j = 1:nc
    S = ones(size(vertices))*diag(ctrxyz(j,:))-1; % Shift matrix
    verts_j = vertices+S;
    hcube(j) = patch('Vertices',verts_j,'Faces',faces,...
        'Parent',hax,...
        'FaceColor',faceClr(j,:),...
        'EdgeColor',edgeClr(j,:),...
        'FaceAlpha',1,'EdgeAlpha',1,...
        'UserData', [ctrxyz(j,1:2),objid(j)]); % [slice, phase, object]
    %             if size(clr,1) == 1
    %                 set(hcube(j),'FaceColor',clr)
    %             else
    %                 set(hcube(j),'FaceColor',clr(j,:),'EdgeColor',clr(j,:)/1.5)
    %             end
end


% ------------------------------------------------------------------------
function cubeBDF(hCube,eventdata)
% The selected hCube may not be the intuitive selection because of Matlab's
% stacking order.  So we'll work out for ourselves which cube was clicked: 

spo = get(hCube,'UserData'); % [slice, phase, object]

handles = guidata(hCube);
hf = handles.figure1;
hax = handles.axes1;
hcubes = findobj(hax,'Type','Patch');
SPO = cell2mat(get(hcubes,'UserData'));

% Reduce data sets to include only this slice (slices don't seem to cause
% so much of a problem):
sInds = SPO(:,3) == spo(3);
SPO = SPO(sInds,:);
hcubes = hcubes(sInds);

% Find neighbouring cubes (inclusive):
idx = find( (abs(spo(1)-SPO(:,1)) <= 1)  & (abs(spo(2)-SPO(:,2)) <= 1));

% Reduce data set:
SPO = SPO(idx,:);
hcubes = hcubes(idx);

% Now how do we fix up the selection of hCube...?
%keyboard

% lw = get(hCube,'LineWidth');
% set(hCube,'LineWidth',2.5);
% pause(0.5);
% set(hCube,'LineWidth',lw);

tobj = timer('ExecutionMode','singleShot','StartDelay',0.5,...
    'StartFcn',@setline,...
    'TimerFcn',@setline,...
    'StopFcn',@setline,...
    'UserData',hCube);
start(tobj);

% Do a silent update to the display of the main program so that it doesn't
% re-trigger the listener to refresh the navigator:
callerHdls = guidata(handles.caller);
updateSlice(callerHdls,spo(1),spo(2),'silent')


% ------------------------------------------------------------------------
function setline(tobj,event)
% Highlights and un-highlights seleced cube(s)
hCube = get(tobj,'UserData');
handles = guidata(hCube);
switch event.Type
    case 'StartFcn'
        hHi = getappdata(handles.figure1,'Highlighted');
        set(hHi,'LineWidth',0.5)
        set(hCube,'LineWidth',2.5)
        setappdata(handles.figure1,'Highlighted',hCube)
        
    case 'TimerFcn'
        set(hCube,'LineWidth',0.5)
        setappdata(handles.figure1,'Highlighted',[])
        
    case 'StopFcn'
        delete(tobj)
        
end



% % ------------------------------------------------------------------------
% function restack(varargin)
% 
% if nargin == 2
%     [~,eventData] = varargin{:};
%     % Only act if the button is turned off; ie, once the user has finished
%     % rotating:
%     if isequal('on', get(eventData,'NewValue'))
%         return
%     end
%     handles = guidata(eventData.AffectedObject);
% else
%     keyboard
% end
% 
% hCubes = findobj(handles.axes1,'type','patch');
% spo = cell2mat([get(hCubes,'UserData')]);
% 
% % Now we first sort
% X0 = spo;
% I0 = (1:size(X0,1))';
% [X1,I1] = sortrows(X0,1);
% I1 = I0(I1);
% [X2,I2] = sortrows(X1,1);
% I2 = I1(I2);
% [X3,I3] = sortrows(X2,3);
% I3 = I2(I3);
% 
% tic
% hCubes = hCubes(I3);
% for j = 1:numel(hCubes)
%     uistack(hCubes(j),'bottom')
% end
% toc
% keyboard


% ------------------------------------------------------------------------
function selectWBDF(hf,~)

[hObj,phase,slice] = getObjectPhase(hf)

if isempty(hObj)
    % Did not select a valid object
    return
end
handles = guidata(hf);
callerHdls = guidata(handles.caller);

% Update the display of the main program:
updateSlice(callerHdls,slice,phase)


% ------------------------------------------------------------------------
function [hObj,phase,slice] = getObjectPhase(hf)

hObj = [];
phase = [];
slice = [];

% Only act on double click?
%if ~isequal('open',get(hf,'SelectionType'))
%    return
%end

% Using the figure's 'CurrentObject' property doesn't work so well with
% BAR3 objects.  Let's do it a different way
%cobj = get(hf,'CurrentObject');
%get(cobj,'tag')


handles = guidata(hf);

% Only use clicks that are within axes bounds:
if ~clickedInAxes(handles.axes1)
    return
end

% Get the clicked point:
cp = get(handles.axes1,'CurrentPoint');

% Now find the intersection of that clicked point with the planes of all
% the bar objects
hb = handles.hbars;
no = numel(hb);
Is = NaN(no,3);             % Intersection matrix
for j = 1:no
    x = get(handles.hbars(j),'XData');
    x = nanmean(x(:));
    Is(j,:) = plane_line_intersect([1 0 0],[x 0 0],cp(1,:),cp(2,:));
end

% Now we have intersections in order of ascending x-value
% This is fine if we are looking at the front of the plot - we will take
% the first object that it intersects.  But if we are looking at the rear
% of the plot we need to sort in reverse order.  We can tell which side
% we're on by looking at the azimuth value:
%   azimuth = -ve  ==>  veiwing from front
%   azimuth = -ve  ==>  viewing from rear
[az,~] = view(handles.axes1);
if az > 0
    Is = flipud(Is); % Reverse if necessary
    hb = flipud(hb);
end

% Get which phase the point intersected each object at, and
%   how many slices it is expecting at each object
phases = round(Is(:,2));
slices = floor(Is(:,3))+1; % ie, 0-1 is slice #1, 1-2 is slice #2, etc

objId = [];
for j = 1:no
    if ( phases(j) < 1 ) || ( phases(j) > size(handles.NS,1))
        continue
    end
    actNslices = handles.NS(phases(j),j);    % actual number of slices here
    if actNslices >= slices(j) && slices(j) > -0.1
        objId = j;
        break
    end
end

% Did we pick one?
if isempty(objId)
    % No...
    return
end

% Got object:
hObj = hb(objId);
%tag = get(hObj,'tag')
phase = phases(objId);           % Selected phase
slice = round(slices(objId));   % Selected slice



% ------------------------------------------------------------------------
function tf = clickedInAxes(ha)
% Check if the click was inside the axes or outside.  See documentation for
% axes property 'CurrentPoint':
cp = get(ha,'CurrentPoint');
lims = [xlim' ylim' zlim'];
% Avoid rounding errors:
f = 1e10;
cp   = round(cp*f)/f;
lims = round(lims*f)/f;
% Check if CurrentPoint is within the axis box:
tf = all( cp(1,:) >= lims(1,:) ) & all( cp(2,:) <= lims(2,:) );


% ------------------------------------------------------------------------
function [I,check]=plane_line_intersect(n,V0,P0,P1)
%plane_line_intersect computes the intersection of a plane and a segment(or
%a straight line)
% Inputs:
%       n: normal vector of the Plane
%       V0: any point that belongs to the Plane
%       P0: end point 1 of the segment P0P1
%       P1:  end point 2 of the segment P0P1
%
%Outputs:
%      I    is the point of interection
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3=>the intersection lies outside the segment P0P1
%
% Example:
% Determine the intersection of following the plane x+y+z+3=0 with the segment P0P1:
% The plane is represented by the normal vector n=[1 1 1]
% and an arbitrary point that lies on the plane, ex: V0=[1 1 -5]
% The segment is represented by the following two points
% P0=[-5 1 -1]
%P1=[1 2 3]
% [I,check]=plane_line_intersect([1 1 1],[1 1 -5],[-5 1 -1],[1 2 3]);

%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%If you have any comments or face any problems, please feel free to leave
%your comments and i will try to reply to you as fast as possible.

I=[0 0 0];
u = P1-P0;
w = P0 - V0;
D = dot(n,u);
N = -dot(n,w);
check=0;
if abs(D) < 10^-7        % The segment is parallel to plane
    if N == 0           % The segment lies in plane
        check=2;
        return
    else
        check=0;       %no intersection
        return
    end
end

%compute the intersection parameter
sI = N / D;
I = P0+ sI.*u;

if (sI < 0 || sI > 1)
    check= 3;          %The intersection point  lies outside the segment, so there is no intersection
else
    check=1;
end



function manual_image_calibration(ax,pixelspacing)
% Calibrate any normal image using a manual scale
%
% This function returns immediately after creating all the necessary
% graphics elements
%
%
% Usage:
%
% Initialising a manual calibration:
%
%   [...] = manual_image_calibration(ax)
%
%
% Displaying a previous calibration for editing / review
%
%   [...] = manual_image_calibration(ax,pix_per_mm)
%
%
% NOTE:   mm_per_pixels = 1/PixelSpacing

I = get( findobj(ax,'type','image'), 'CData');


if ~exist('pixelspacing','var') || isempty(pixelspacing) || ~isfinite( mean(pixelspacing) )
    pixelspacing = NaN; % Flag to get set to default later
    
else
    % input will be [x;x] but we only want one number, so take mean:
    pixelspacing = mean(pixelspacing);
        
end

[xe,ye,pix_per_inc] = configure_line(I);

N = 10;
x = linspace( xe(1),xe(2), N );
y = linspace( ye(1),ye(2), N );

plot_static_line(ax,x,y);

show_dialog(ax,N, pixelspacing, pix_per_inc)


% ------------------------------------------------------------------------
function show_dialog(ax, N, pixelspacing, pix_per_inc)

% Calculation of parameters & configuring defaults

hf = ancestor(ax,'Figure');

% Dialog setup:
fpsn = get(hf,'OuterPosition');
w = 200;
h = 260;
psn(1) = fpsn(1)+fpsn(3)+2;
psn(2) = fpsn(2)+fpsn(4) - h;
psn(3:4) = [w h];


fig = figure('Toolbar','none',...
    'MenuBar','none',...
    'Units','pixels',...
    'NumberTitle','off',...
    'HandleVisibility','off',...
    'Name','Calibration',...
    'OuterPosition', psn,...
    'CloseRequestFcn',@fig_closereq);

setappdata(fig,'target_axes',ax);

fs = 11;

y = h-60;
hcalLine = uicontrol(fig,'Style','RadioButton',...
    'String','Calibration Line',...
    'Tag','calLine_radiobutton',...
    'FontSize',fs,...
    'Position',[20 y w-40 20],...
    'Callback',@radio_callback);

y = y-25;
hpts = uicontrol(fig,'Style','Edit',...
    'String',num2str(N),...
    'Tag','npoints_edit',...
    'FontSize',fs,...
    'BackgroundColor',[1 1 1],...
    'Position',[45 y 40 22],...
    'Callback',@edit_points_callback);

uicontrol(fig,'Style','Text',...
    'String','points',...
    'FontSize',fs,...
    'HorizontalAlignment','left',...
    'BackgroundColor',get(fig,'Color'),...
    'Position',[87 y 80 18]);

y = y-25;
hspc = uicontrol(fig,'Style','Edit',...
    'String',sprintf('%2.1f',1),...
    'Tag','spacing_edit',...
    'FontSize',fs,...
    'BackgroundColor',[1 1 1],...
    'Position',[45 y 40 22],...
    'Callback',@edit_spacing_callback);

uicontrol(fig,'Style','Text',...
    'String','mm / increment',...
    'FontSize',fs,...
    'HorizontalAlignment','left',...
    'BackgroundColor',get(fig,'Color'),...
    'Position',[87 y 80 18]);

y = y-40;
hdirect = uicontrol(fig,'Style','RadioButton',...
    'String','Direct',...
    'Tag','direct_radiobutton',...
    'FontSize',fs,...
    'Position',[20 y w-40 20],...
    'Callback',@radio_callback);

y = y-25;
hpixspac = uicontrol(fig,'Style','Edit',...
    'String',num2str(pixelspacing),...
    'Tag','pixelspacing_edit',...
    'FontSize',fs,...
    'BackgroundColor',[1 1 1],...
    'Position',[45 y 60 22],...
    'Callback',@edit_pixelspacing_callback);
uicontrol(fig,'Style','Text',...
    'String','"PixelSpacing"',...
    'FontSize',fs,...
    'HorizontalAlignment','left',...
    'BackgroundColor',get(fig,'Color'),...
    'Position',[107 y 80 18]);

y = y-25;
uicontrol(fig,'Style','Edit',...
    'String',num2str(1/pixelspacing),...
    'Tag','pixpermm_edit',...
    'FontSize',fs,...
    'BackgroundColor',[1 1 1],...
    'Position',[45 y 60 22],...
    'Callback',@edit_pixpermm_callback);
uicontrol(fig,'Style','Text',...
    'String','pixels/mm',...
    'FontSize',fs,...
    'HorizontalAlignment','left',...
    'BackgroundColor',get(fig,'Color'),...
    'Position',[107 y 80 18]);

y = y-50;
uicontrol(fig,'Style','Pushbutton',...
    'String','Save Calibration',...
    'Tag','calibrate_pushbutton',...
    'Callback',@save_calibration,...
    'FontSize',fs,...
    'FontWeight','bold',...
    'Position',[65 y 105 30])

disp('stuff here not fully resolved')


% If pixelspacing is not specified (NaN), then select "Calibration Line" as
% default and update pixelspacing from that:
if isnan(pixelspacing)
    target = hcalLine;
    update_pixel_spacing_text(hspc)
    
else
    target = hdirect;
    
end
%set(target,'Value',1)
radio_callback(target,[])


% Add listener to update pixelspacing texts:
ax = getappdata(fig,'target_axes');
he = findobj(ax,'Tag','endpoints');
addlistener( he,'YData', 'PostSet',@(a,b)line_listener(a,b,fig) );


% ------------------------------------------------------------------------
function line_listener(~,~,cal_figure)
% This listener responds to movement of the line and updates pixel spacing
% texts accordingly:
hObject = findobj( get(cal_figure,'Children'), 'Tag', 'spacing_edit' );
update_pixel_spacing_text(hObject)


% ------------------------------------------------------------------------
function tf = equal_with_tol(a,b,tol)
tf = abs(a-b) < tol;


% ------------------------------------------------------------------------
function radio_callback(hObject,~)

% Mutual exclusivitiy:
fig = ancestor(hObject,'Figure');
kids = get(fig,'Children');
hrad = findobj(kids,'-depth',0,'-regexp','tag','_radiobutton');

set(hrad(hrad ~= hObject),'Value',0)
set(hObject,'Value',1)

% Hide / show calibration line

% Enable/disable edit boxes
line_edits = findobj(kids,'Tag','npoints_edit','-or','Tag','spacing_edit');
direct_edits = findobj(kids,'Tag','pixelspacing_edit','-or','Tag','pixpermm_edit');

line = get_static_line(fig);

switch get(hObject,'Tag')
    case 'calLine_radiobutton'
        set(line,'Visible','on')
        set(line_edits,  'Enable','on')
        set(direct_edits,'Enable','off')
    case 'direct_radiobutton'
        set(line,'Visible','off')
        set(line_edits,  'Enable','off')
        set(direct_edits,'Enable','on')
end

% ------------------------------------------------------------------------
function pixelspacing = cal_line_pixelspacing(fig)
% Calculate pixel spacing from calibration line:

kids = get(fig,'Children');

% Find the line's endpoints:
ax = getappdata(fig,'target_axes');
he = findobj(ax,'Tag','endpoints');
xe = get(he,'XData');
ye = get(he,'YData');

% Number of points on the line:
n = str2double( get( findobj(kids,'tag','npoints_edit'), 'String') );

% Lenght of an increment:
spacing = str2double( get( findobj(kids,'tag','spacing_edit'), 'String') );

% Pixel length of the line:
len_pix = sqrt( (xe(2)-xe(1)).^2 + (ye(2)-ye(1)).^2 );

% Real length of the line
len_mm = spacing * (n-1);

% pixelspacing = mm / pix
pixelspacing = len_mm / len_pix;

% ------------------------------------------------------------------------
function edit_pixpermm_callback(hObject,~)

update_pixel_spacing_text( hObject )


% ------------------------------------------------------------------------
function edit_pixelspacing_callback(hObject,~)

update_pixel_spacing_text( hObject )


% ------------------------------------------------------------------------
function edit_points_callback(hObject,~)
n = str2num( get(hObject,'String') ); %#ok<ST2NM>
if ~isfinite(n)
   n = 10;
end
n = validate_in_range(n,2,inf,true);
set(hObject,'String', num2str(n))

ax = getappdata(gcbf,'target_axes');
update_static_line(ax,n)

update_pixel_spacing_text( hObject );


% ------------------------------------------------------------------------
function edit_spacing_callback(hObject,~)
v = str2num( get(hObject,'String') ); %#ok<ST2NM>
if ~isfinite(v)
    v = 1;
end
v = validate_in_range(v,eps,inf,false);
set(hObject,'String', num2str(v) );

update_pixel_spacing_text( hObject );


% ------------------------------------------------------------------------
function update_pixel_spacing_text(hObject)
% Update all displays of pixelspacing or derivatives:


switch get(hObject,'Tag')
    case 'pixelspacing_edit'
        pixelspacing = str2double( get(hObject,'String') );
        pixpermm = 1/pixelspacing;
        
    case 'pixpermm_edit'
        pixpermm = str2double( get(hObject,'String') );
        pixelspacing = 1/pixpermm;
        
    otherwise
        % Calculate from calibration line
        fig = ancestor( hObject, 'figure' );
        pixelspacing = cal_line_pixelspacing( fig );
        pixpermm = 1/pixelspacing;
end

hf = ancestor(hObject,'figure');
kids = get(hf,'Children');

set( findobj(kids,'Tag','pixelspacing_edit'), 'String', num2str( pixelspacing ) );
set( findobj(kids,'Tag','pixpermm_edit'),     'String', num2str( pixpermm ) );


% ------------------------------------------------------------------------
function val = validate_in_range(val,minval,maxval,integer)
if integer
    val = round(val);
end
val = max([val minval]);
val = min([val maxval]);



% ------------------------------------------------------------------------
function [hg,hl,he] = plot_static_line(ax,x,y)

xe = x([1,end]);
ye = y([1,end]);

hold(ax,'on')
hg = hggroup('Parent',ax,'Tag','calibration_line','HitTest','off');
hl = plot(ax,x,y,'g+','LineStyle','-','Tag','discretised_line','Parent',hg);
he = plot(ax,xe,ye,'r+','MarkerSize',14,'Tag','endpoints','Parent',hg);

% Configure interactivity:
set(hl,'ButtonDownFcn',@button_down)
set(he,'ButtonDownFcn',@button_down)


%{
% Add responsive cursors if possible:
if exist('iptSetPointerBehavior','file') == 2
    hf = get(ax,'Parent');
    enterFcn = @(figHandle, currentPoint)setpointer(figHandle,'fleur');
    iptSetPointerBehavior([hl;he],enterFcn);
    iptPointerManager( hf );
end
%}



% ------------------------------------------------------------------------
function save_calibration(hObject,eventdata)
% This function saves the calibration to the ImageStack object in the
% "Images" property of handles, then destroys the interactive tools.
%

% Get graphics handles:
fig = ancestor(hObject,'figure');

ax = getappdata(fig,'target_axes');
hf = get(ax,'Parent');

pixelspacing = str2double( get( findobj(get(fig,'children'),'Tag','pixelspacing_edit'), 'String') );

handles = guidata(hf);

% Push data back to parent GUI:
handles.Images.s = [1;1]*pixelspacing;
guidata(hf,handles)

% Destroy figure:
fig_closereq(fig)





% ------------------------------------------------------------------------
function hg = get_static_line(fig)
ax = getappdata(fig,'target_axes');
hg = findobj(ax,'Type','hggroup','Tag','calibration_line');

% ------------------------------------------------------------------------
function update_static_line(ax,N)
he = findobj(ax,'Tag','endpoints');
xe = get(he,'XData');
ye = get(he,'YData');

xl = linspace(xe(1),xe(end),N);
yl = linspace(ye(1),ye(end),N);

hl = findobj(ax,'Tag','discretised_line');
set(hl,'XData',xl,'YData',yl)


% ------------------------------------------------------------------------
function fig_closereq(fig,~)

% Remove the calibration line:
try %#ok<TRYNC>
    hg = get_static_line(fig);
    delete(hg);
end

closereq


% ------------------------------------------------------------------------
%       Interactive callbacks
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
function button_down(obj,~)

fig = gcf;
ax = gca;

% Add the listeners
add_listener(fig,'WindowButtonMotionEvent',@motion_fcn)
add_listener(fig,'WindowButtonUpEvent',@button_up)

% Store cursor position and selected object:
cp = get(ax,'CurrentPoint');
setappdata(ax,'CursorLocation',cp(1,1:2));
setappdata(ax,'SelectedObject',obj);


% ------------------------------------------------------------------------
function button_up(fig,~)

% Remove the listeners
remove_listener(fig,'WindowButtonMotionEvent',@motion_fcn)
remove_listener(fig,'WindowButtonUpEvent',@button_up)


% ------------------------------------------------------------------------
function motion_fcn(fig,eventData)

% Here we re-draw the line using the same number of points, but by moving
% the end-points & re-creating & re-setting the x & y data

% --- Persistents for speed --- %
persistent ax hl he

valid_handle = @(h) ~isempty(h) && ishandle(h);

if ~valid_handle(ax)
    ax = findobj( get(fig,'Children'),'-depth',0,'Type','axes','Tag','axes1');
end
if ~valid_handle(hl)
    hl = findobj(ax,'Type','line','Tag','discretised_line');
end
if ~valid_handle(he)
    he = findobj(ax,'Type','line','Tag','endpoints');
end
% --- ..................... --- %
 

xe = get(he,'XData');
ye = get(he,'YData');


% When listenting to a figure's WindowButtonMotionEvent, we do not get
% access to an updated cursor point via
%   get(ax,'CurrentPoint')
% Instead, we have to use the current point from eventdata:
currentPointFig = eventData.CurrentPoint;
% However, this current point comes to us in figure coordinates, and in
% pixels.  So we need to convert to axis coordinates, in data units
currentPoint = figPostion2axisDataPosition(fig,ax,currentPointFig);

xp = currentPoint(1);
yp = currentPoint(2);

% Now switch the activity based on which object is being manipulated:
%        line => pan entire hggroup
%   endpoints => move selected endpoints & re-calc intermediate points
obj = getappdata(ax,'SelectedObject');
if obj == hl
    dxy = currentPoint - getappdata(ax,'CursorLocation');
    setappdata(ax,'CursorLocation',currentPoint);   
    
    xe = xe + dxy(1);
    ye = ye + dxy(2);
    
    xl = get(hl,'XData') + dxy(1);
    yl = get(hl,'YData') + dxy(2);
    
elseif obj == he
    
    d = sqrt((xe-xp).^2 + (ye-yp).^2);  % [1x2] distances to the two endpoints
    ptid = d == min(d);                 % point under cursor
    
    xe(ptid) = currentPoint(1);
    ye(ptid) = currentPoint(2);
    
    % Now re-calculate the line
    N = numel(get(hl,'XData'));
    xl = linspace(xe(1),xe(end),N);
    yl = linspace(ye(1),ye(end),N);
    
else
    warning('Unhandled object') %#ok<WNTAG>
    return
end

% Now update the data:
set(he,'XData',xe,'YData',ye)
set(hl,'XData',xl,'YData',yl)


% ------------------------------------------------------------------------
function pt = figPostion2axisDataPosition(fig,ax,psn)
% PSN is [1-by-2] location within FIG, in pixels
%
% AX is the axis in which we wish to get the location.

% Convert to normalized axis coordinates:
AxisPos = hgconvertunits(fig,get(ax,'position'),get(ax,'units'), ...
    'pixels',fig);
ptAxis = (psn - AxisPos(1:2)) ./ AxisPos(3:4);

NormX = ptAxis(1);
NormY = ptAxis(2);

% Check for inverted axes:
isRevX = strcmp(get(ax,'XDir'),'reverse');
isRevY = strcmp(get(ax,'YDir'),'reverse');
if isRevX
    NormX = 1-NormX;
end
if isRevY
    NormY = 1-NormY;
end

Xlim = get(ax,'Xlim');      % X limits
Ylim = get(ax,'Ylim');      % Y limits

pt = [NormX NormY].* [Xlim(2)-Xlim(1) Ylim(2)-Ylim(1)] + [Xlim(1) Ylim(1)];



% ------------------------------------------------------------------------
function add_listener(obj,varargin)

listeners = getappdata(obj,'listeners');

% Instantiate listener
li = handle.listener(obj,varargin{:});

% Append, handling the empty case
if isempty(listeners)
    listeners = li;
else
    listeners(end+1) = li;
end

% Store for later:
setappdata(obj,'listeners',listeners)



% ------------------------------------------------------------------------
function remove_listener(obj,eventType,fun)

li = getappdata(obj,'listeners');

event_eq = strcmp( get(li,'EventType'), eventType );
for j = numel(li) : -1 : 1
    cbk_eq(j,1) = isequal(li(j).Callback, fun); 
end

sel = event_eq & cbk_eq;

delete(li(sel))

li(sel) = [];

setappdata(obj,'listeners',li);


% ------------------------------------------------------------------------
function [xe,ye,ppi] = configure_line(I)
% Create the data & scale for the line that will be displayed
%
%  xe - xdata of endpoints of the line
%  ye - ydata of endpoints of the line
% ppi - pixels per increment 

r = size(I,1);
c = size(I,2);

ye = [1 1] * round(2/3*r);
xe = round( [1/6 5/6]*(c-1) + 1 );

ppi = NaN;


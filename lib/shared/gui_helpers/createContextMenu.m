function hcmenu = createContextMenu(opt)
%CREATECONTEXTMENU Create contextual menus for different items
%
% OPT is a menu type:
%   'cloud'
%   'trace'

hcmenu = uicontextmenu; % Menu top level

switch lower(opt)

    case 'trace'
        
        set(hcmenu,'Tag','Trace context menu')
        
        % Colour Sub-menu
        colourMenu(hcmenu);
        
        % Delete item
        uimenu(hcmenu,...
            'Label',    'Delete ROI',...
            'Separator','on',...
            'Callback', @deleteROI);
        
    case {'cloud','staticcloud'}
        set(hcmenu,'Tag','Cloud context menu')
        % Create a context menu for point clouds:
        % Provides:
        %   - Colour options
        %   - Delete 
        
        % Graphics: 
        colourMenu(hcmenu);     % Colour options
        lineStyleMenu(hcmenu);  % Line style options
        lineWidthMenu(hcmenu);  % Line width options
        markerMenu(hcmenu);     % Marker style options
        markerSizeMenu(hcmenu); % Marker size options
        
        % Manually Adjust Position
        %uimenu(hcmenu,'Label','Adjust Position','Callback',@adjustCloudPsn);
        
        % Optimise Position
        %uimenu(hcmenu,'Label','Optimise Position','Callback',@optimiseCloudPsn);
        
        % Hide item:
        %uimenu(hcmenu,'Label','Hide','Separator','on','Callback',@hideCloud);
        
        % Delete item:
        %uimenu(hcmenu,'Label','Delete','Separator','on','Callback',@deleteCloud);
        
    otherwise
        warning(['No menu option for specifier "' opt '".']);
        
end




% ------------------------------------------------------------------------
function item1 = colourMenu(hparent)
item1 = uimenu(hparent,'Label','Colour','Callback',@updateCheck);
cspec = ColourSpec();       % Get the ColourSpec map
cspec = cspec(:,1);            % Just keep the long names
for j = 1:numel(cspec)      % Loop to build all colour entries
    uimenu(item1,'Label',cspec{j},'Callback',@setColour);
end


% ------------------------------------------------------------------------
function item1 = lineStyleMenu(hparent)
item1 = uimenu(hparent,'Label','Line Style','Callback',@updateCheck);
lspec = LineSpec();     % Get the LineSpec map
lspec(:,2) = [];        % Keep only the names
for j = 1:numel(lspec)
    uimenu(item1,'Label',lspec{j},'Callback',@setLineStyle);
end


% ------------------------------------------------------------------------
function item1 = lineWidthMenu(hparent)
item1 = uimenu(hparent,'Label','Line Width','Callback',@updateCheck);
opts = {'0.5','1.0','2.0','3.0','4.0','5.0','6.0',...
    '7.0','8.0','9.0','10.0','11.0','12.0'};
for j = 1:numel(opts)
    uimenu(item1,'Label',opts{j},'Callback',@setLineWidth)
end
    

% ------------------------------------------------------------------------
function item1 = markerMenu(hparent)
item1 = uimenu(hparent,'Label','Marker','Callback',@updateCheck);
opts = {'+','o','*','.','x','square','diamond','v','^','>','<',...
    'pentagram','hexagram','none'};
for j = 1:numel(opts)
    uimenu(item1,'Label',opts{j},'Callback',@setMarker);
end


% ------------------------------------------------------------------------
function item1 = markerSizeMenu(hparent)
item1 = uimenu(hparent,'Label','Marker Size','Callback',@updateCheck);
opts = {'2','4','6','8','10','12','18','24','48'};
for j = 1:numel(opts)
    uimenu(item1,'Label',opts{j},'Callback',@setMarkerSize);
end


% ------------------------------------------------------------------------
function updateCheck(hmenu,~)
%UPDATECHECK Update the check marks on submenu

h = gco;    % current object

menuName = get(hmenu,'Label');
items = get(hmenu,'Children');  % ITEMS will always list from bottom up.
items = flipud(items);          % reverse for convenience
lbls  = get(items,'Label');

switch lower(menuName)
    case 'colour'
        clist = ColourSpec('torgb',lbls);
        tf = all( clist == repmat(get(h,'Color'),size(clist,1),1), 2);
        
    case 'line style'
        map = LineSpec();
        tf = strcmp(get(h,'LineStyle'),map(:,2));
        
    case 'line width'
        tf = get(h,'LineWidth') == str2double(lbls);
        
    case 'marker'
        tf = strcmpi(get(h,'Marker'),lbls);
        
    case 'marker size'
        tf = get(h,'MarkerSize') == str2double(lbls);
        
    otherwise
        tf = false(size(lbls));
        
end
set(items,'Checked','off');
set(items(tf),'Checked','on');


% ------------------------------------------------------------------------
function setLineStyle(hmenu,~)
%SETLINESTYLE Set line style of current object
%
h = gco;    % Line handle
lineStr = get(hmenu,'Label');
sty = LineSpec({lineStr});
set(h,'LineStyle',sty{:});


% ------------------------------------------------------------------------
function setLineWidth(hmenu,~)
%SETLINEWIDTH Set line width of line object
h = gco;
lw = str2double(get(hmenu,'Label'));
set(h,'LineWidth',lw)


% ------------------------------------------------------------------------
function setColour(hmenu,~)
%SETCOLOR Set colour of current object.

% The other option would be to call UISETCOLOR

h = gco;	% Line handle

clrString = get(hmenu,'Label');

% Update the graphics object:
set(h,'Color',clrString);

% Get handles to do other stuff:
handles = guidata(hmenu);

% If it's a ROI tracing, then update the data structure:
if isfield(handles,'traces')
    tf = strcmpi(get(h,'Tag'),{handles.traces.Tag});
    if ~isempty(tf) && any(tf)
        handles.traces(tf).Color = clrString;
        guidata(hmenu,handles)
    end
end

% If it's a static cloud, then update the current axes ColorOrder
if isfield(handles,'Models')
    tf = strcmpi(get(h,'Tag'),{handles.Models.Tag});
    if ~isempty(tf) && any(tf)
        cOrder = get(gca,'ColorOrder'); % Get current ColorOrder list
        cOrder(tf,:) = get(h,'Color');  % Edit the tibia's one to match user selection
        set(gca,'ColorOrder',cOrder)    % Update axes property
    end
end
    

% ------------------------------------------------------------------------
function setMarker(hmenu,~)
%SETMARKER Set marker of current object
%
h = gco;    % Line object handle
set(h,'Marker',get(hmenu,'Label'));


% ------------------------------------------------------------------------
function setMarkerSize(hmenu,~)
%SETMARKER Set marker of current object
%
h = gco;    % Line object handle
set(h,'MarkerSize',str2double(get(hmenu,'Label')));


% ------------------------------------------------------------------------
function deleteROI(hmenu,~)
handles = guidata(hmenu);

ht = gco;   % Trace handle

% Remove selected trace:
tf = strcmpi(get(ht,'Tag'),{handles.traces.Tag});
handles.traces(tf) = [];

% Update handles - this must come before updateSlices() because
% that function trigerrs other functions which are dependent on
% handles being up to date.
guidata(hmenu,handles)

% Redraw traces:
updateSlice(handles)


% ------------------------------------------------------------------------
function deleteCloud(hObject,~)
% DELETECLOUD Contextual menu callback to delete point clouds
%
% hObject is the handle to the current uimenu item

handles = guidata(hObject);
hcld = gco;

% Get tags of all clouds & associated CS:
hallCld = [handles.HiResClouds.hObj];
hallCS  = [handles.HiResClouds.hObjCS];

% Delete graphics objects:
cid = find(hallCld==hcld);
delete(hcld)
delete(hallCS(cid))

% Remove stored data:
handles.HiResClouds(cid) = [];

% Update handes:
guidata(hObject,handles);


% ------------------------------------------------------------------------
function hideCloud(hObject,~)
% HIDECLOUD Contextual menu callback to temporarily remove point cloud from
% current display
%
% hObject is the handle to the current uimenu item

hcld = gco;
set(hcld,'Visible','off')


function varargout = ModelLoader(varargin)
% MODELLOADER MATLAB code for ModelLoader.fig
%      MODELLOADER, by itself, creates a new MODELLOADER or raises the existing
%      singleton*.
%
%      H = MODELLOADER returns the handle to a new MODELLOADER or the handle to
%      the existing singleton*.
%
%      MODELLOADER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODELLOADER.M with the given input arguments.
%
%      MODELLOADER('Property','Value',...) creates a new MODELLOADER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ModelLoader_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ModelLoader_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ModelLoader

% Last Modified by GUIDE v2.5 18-Oct-2012 14:49:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ModelLoader_OpeningFcn, ...
    'gui_OutputFcn',  @ModelLoader_OutputFcn, ...
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
end


% ------------------------------------------------------------------------
function ModelLoader_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before ModelLoader is made visible.
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ModelLoader (see VARARGIN)

% Choose default command line output for ModelLoader
handles.output = hObject;

% Remove all children & build from scratch:
delete(get(handles.figure1,'Children'))

% Sticky path:
handles.userPath = [pwd filesep];

% Position gui over the top of calling gui
positionOver(handles.figure1,gcbf);

% Make UITree - this also updates handles:
buildGui(handles,varargin{:});

% UIWAIT makes ModelLoader wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end


% ------------------------------------------------------------------------
function varargout = ModelLoader_OutputFcn(hObject, eventdata, handles)
% --- Outputs from this function are returned to the command line.
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end


% ------------------------------------------------------------------------
function buildGui(handles,models)
%BUILDGUI Main function which sets up the gui and builds the uitree by
%populating with input data

% --------- Set up buttons:
w = 0.3;    % width
h = 0.05;   % height
sep = 0.02; % vertical separation

basic = [1-w-sep, 1-h-sep/2, w, h];
handles.AddButton = uicontrol('String','Add Object',...
    'Enable','on',...
    'Units','normalized',...
    'Position',basic - [0, sep, 0, 0],...
    'Callback',@AddButton_Cbk);
handles.RemoveButton = uicontrol('String','Remove Object',...
    'Enable','off',...
    'Units','normalized',...
    'Position',basic - [0, 2*sep+h, 0, 0],...
    'Callback',@RemoveButton_Cbk);

handles.SmartAddButton = uicontrol('String','Smart Add Files',...
    'Enable','on',...
    'Units','normalized',...
    'Position', basic - [0 4*sep+2*h, 0, 0],...
    'Callback',@SmartAddButton_Cbk);

basic = [1-w-sep, 0+sep/2, w, h];
handles.OkButton = uicontrol('String','Ok',...
    'Units','normalized',...
    'Position',basic + [0, sep, 0, 0],...
    'Callback',@OkButton_Cbk);
handles.CancelButton = uicontrol('String','Cancel',...
    'Units','normalized',...
    'Position',basic + [0 2*sep+h, 0, 0],...
    'Callback',@CancelButton_Cbk);

dw = 0.05;
handles.ClearButton = uicontrol('String','',...
    'Enable','on',...
    'Units','normalized',...
    'Position', basic + [dw/2, 3*sep+4*h, -dw, 0],...
    'Callback',@ClearButton_Cbk);



% --------- Build tree:
import javax.swing.*
import javax.swing.tree.*;

% Create java icons:
javaImage_checked = getIcon('checked');
javaImage_unchecked = getIcon('unchecked');

 
% Images javaImage_checked/unchecked are assumed to have the same width
iconWidth = javaImage_unchecked.getWidth;

% Create top node
warning('off','MATLAB:uitreenode:DeprecatedFunction')
rootNode = uitreenode('v0','root', 'Models', [], 0);

% Set treeModel
treeModel = DefaultTreeModel( rootNode );

% Create the tree
tree = uitree('v0');
tree.setModel( treeModel );

% Store the figure handle in the UserObject of the root node.  This way we
% can get access to handles, and then to everything else
rootNode.setUserObject( handles.figure1 );

% Some things supported by jtree, not the uitree:
jtree = handle(tree.getTree,'CallbackProperties');

% Store these:
handles.tree = tree;
handles.jtree = jtree;

% Create children if necessary, and update guidata
if nargin > 1 && ~isempty(models)
    pset = {models(1).HiRes(:).Path, ... % by concatenating these
            models(1).LoRes(:).Path};    % we'll ensure we get a path
    pth = pset{1};
    handles.userPath = pth(1:find(pth==filesep,1,'last'));
    
    % Populate the tree:
    populateTree(models);
end

% Now update guidata.  This covers 2 cases:
%   - Nargin > 1 case (above):  handles has been modified by populateTree
%   - Nargin = 0 case:          handles un-modified
guidata(handles.figure1,handles)


% Layout 
%   The figure must be made visible before we can set the size of the tree
%   pane.  Ideally we would prefer this was not so.  The next best thing is
%   that we at least don't want to see the tree pane re-size, so we make in
%   invisible until it has been sized, then make it visible:
set(tree,'Visible',false)
set(handles.figure1, 'Visible','on'), drawnow update
set(tree, 'Units', 'normalized', 'position', [0, 0, 1-w-2*sep, 1],'Visible',false);
set(tree, 'NodeSelectedCallback', @selected_cb );

% Now expand to the level of the categories:
catNodes = cell2mat( findnode(rootNode,'getLevel',2) );
refreshNodeSpace( catNodes );

% Make root the initially selected node:
tree.setSelectedNode( rootNode );

% Enable multi-selection:
tree.setMultipleSelectionEnabled( true );

% Now display the tree
set(tree,'Visible',true)

% Configure buttons:
btnStates(handles,tree,rootNode.getLevel());
 
% MousePressedCallback is not supported by the uitree, but by jtree
set(jtree, 'MousePressedCallback', @mousePressedCallback);
set(jtree, 'KeyPressedCallback',   @keyPressedCallback);
set(jtree, 'ValueChangedCallback', @updateClearButtonText)

depth = []; % variable which works as a persistent variable for updateClearButtonText()
updateClearButtonText()

% Done; now just include helper functions

% --------- Sub-functions & callbacks ------------------------------------
% Note that in these subfunctions, if ever HANDLES is required, it will not
% be up to date unless it is retrieved with guidata()

    % ----------------------------
    function updateClearButtonText(varargin)
        % Update the text for uicontrol ClearButton
        % This button has a two-stage function: if files exist, it will
        % clear all files.  If no files exist it will clear all objects.
        if isequal(depth,rootNode.depth)
            % No change to tree depth -> nothing to do
           return
        end
        depth = rootNode.getDepth;
        if depth == 3
            str = 'Clear all files';
            ev  = 'on';
        elseif depth == 2
            str = 'Clear all objects';
            ev  = 'on';
        else
            str = 'Clear...';
            ev  = 'off';
        end
        set(handles.ClearButton,'String',str,'Enable',ev,'Visible',ev)
        
    end %updateClearButtonText()

    % ----------------------------
    function mousePressedCallback(hTree, eventData) %,additionalVar)
        % if eventData.isMetaDown % right-click is like a Meta-button
        % if eventData.getClickCount==2 % how to detect double clicks
        
        % Get the clicked node
        clickX = eventData.getX;
        clickY = eventData.getY;
        treePath = jtree.getPathForLocation(clickX, clickY);
        % Check if a node was clicked
        if ~isempty(treePath)
            
            % Change the action of the buttons depending on the node level:
            node = treePath.getLastPathComponent;
                        
            % Manage icon changes:
            switch node.getLevel
                case 0  % Root node clicked
                    
                case 1  % Object node clicked
                    
                case 2  % "Static" or "Dynamic" node clicked
                    
                    % Check if the Icon was clicked:
                    if clickX <= (jtree.getPathBounds(treePath).x+iconWidth)
                        nodeValue = node.getValue;
                        % as the value field is the selected/unselected flag,
                        % we can also use it to only act on nodes with these values
                        switch nodeValue
                            case 'selected'
                                setNodeVis(node,'unselected');
                            case 'unselected'
                                setNodeVis(node,'selected');
                        end
                    end
                    
                case 3  % File node clicked
                                        
            end
                    
            
        end
    end %mousePressedCallback()

    % ----------------------------
    function setNodeVis(node,state)
        switch state
            case {'selected', true, 1}
                node.setValue('selected')
                node.setIcon(javaImage_checked)
                jtree.treeDidChange();
            case {'unselected', false, 0}
                node.setValue('unselected')
                node.setIcon(javaImage_unchecked)
                jtree.treeDidChange();
        end
    end %setNodeVis()

    % ----------------------------
    function selected_cb( tree, ~ )
        % The call to getSelectedNodes might occasionally return empty,
        % probably if the tree is still re-drawing. A drawnow call seems to
        % avoid this:  
        drawnow
        node = first(tree.getSelectedNodes);
        if isempty(node)
            lv = -1;
        else
            lv = node.getLevel;
        end
        btnStates(handles,tree,lv) % Set button states

    end %selected_cb()

    % ----------------------------
    function keyPressedCallback(~, eventData)
        keyCode = eventData.getKeyCode;
        switch keyCode
            case {10, 65} % 'return' or 'a'
                % Return key: Rename object
                node = first(tree.getSelectedNodes);
                if isempty(node)
                    return
                elseif (node.getLevel == 1) && (keyCode == 10)
                    % Objects get a rename action:
                    
                    name0 = char( node.getName );
                    
                    dlgTitle = 'Edit Object';
                    relocateDlg(handles.figure1,dlgTitle);
                    name = inputdlg2('Name:',dlgTitle,1,{name0});
                    if isempty( name )
                        return
                    else
                        name = strtrim(name{1});
                    end
                    if isequal(name,name0) || nameInUse( tree, name )
                        return
                    else
                        node.setName( name )
                        refreshNodeSpace( node )
                    end
                    
                elseif any( node.getLevel == [2 3])
                    % Categories & files get an "Add File" action
                    AddButton_Cbk(handles.AddButton,[])
                    
                end
                
        end %switch
        
    end %keyPressedCallback()

    % ----------------------------
    function populateTree(mdls)
        for mj = 1:numel(mdls)  % Loop through models
            oNode = newObject(tree,rootNode,mdls(mj).Tag);
            nc = oNode.getChildCount;
            for cj = 1:nc       % Loop Through Categories: Static | Dynamic
                cNode = oNode.getChildAt(cj-1); % Category node (Static/Dynamic)
                vis = 1;        % In case any categories remain empty
                switch getCatName(cNode)        % Text: 'Static' | 'Dynamic'
                    case 'Static'
                        % Create the node, if cloud exisits:
                        if ~isempty(mdls(mj).HiRes)
                            handles = newFileNode(handles,treeModel,...
                                cNode, 0, mdls(mj).HiRes, mdls(mj).q, mdls(mj).x );
                            vis = mdls(mj).HiRes.Visible;
                        end
                    case 'Dynamic'
                        if ~isempty(mdls(mj).LoRes)
                            for dj = 1:numel(mdls(mj).LoRes)
                                handles = newFileNode(handles,treeModel,...
                                    cNode, dj-1, mdls(mj).LoRes(dj) );
                            end
                            vis = mdls(mj).LoRes(dj).Visible;
                        end
                end %switch
                setNodeVis(cNode,vis)
                updateCatName(tree,cNode)
            end %for
        end %for
        
    end %populateTree()

end %buildGui()


% ------------------------------------------------------------------------
function tf = nameInUse(tree,name)
%NAMEINUSE Check to see if the model name NAME is already in use in TREE
tf = false;

root = tree.getModel().getRoot();   % Root node
objnames = getChildNames(root);

% Test for conflict:
if any(strcmpi(name,objnames))
    dlgTitle = 'Name in use';
    relocateDlg(root.getUserObject,dlgTitle);
    warndlg('This name is already in use. Please use another name',...
        dlgTitle,'modal');
    tf = true;
end

end %nameInUse()


% ------------------------------------------------------------------------
function names = getChildNames(node)
% Return a cell array containing the names all the immediate children.
nobj = node.getChildCount();
% Traul through all level 1 nodes (the object/model nodes):
names = {};
for j = nobj:-1:1
    names{j,1} = char(node.getChildAt(j-1).getName);
end
end %getChildNames()


% ------------------------------------------------------------------------
function relocateDlg(hparent,dtitle)
% This is a complicated mess...
li = [];
hdlg = [];
t = timer('TimerFcn',@reposition,...
    'Period',0.01,...
    'ExecutionMode','fixedSpacing',...
    'TasksToExecute',Inf,...
    'StartDelay',0);
start(t)
    %----------------------------
    function reposition(varargin)
        if isa(varargin{1},'timer')
            set(0,'ShowHiddenHandles','on')
            hdlg = findobj(0,'Type','figure','Name',dtitle);
            if isempty(hdlg)
                return
            end
            li = addlistener(hdlg,'Visible','PostSet',@reposition);
            stop(varargin{1})
            delete(varargin{1})
        else
            delete(li)
            % Need to get this again, because it got changed:
            hdlg = findobj(0,'Type','figure','Name',dtitle);
            positionOver(hdlg,hparent)
            set(0,'ShowHiddenHandles','off')
        end
    end %reposition();
end %relocateDlg()


% ------------------------------------------------------------------------
function btnStates(handles,tree,level)
% Set the String and Enable states of the action buttons

SA.Enable = 'off';
switch level
    case -1 % No selection
        set(handles.AddButton,'Enable','off','String','Add ...')
        set(handles.RemoveButton,'Enable','off','String','Remove ...')
        
    case 0  % Root node selected
        set(handles.AddButton,'Enable','on','String','Add Object')
        set(handles.RemoveButton,'Enable','off','String','Remove Object')
        if tree.getModel.getRoot.getChildCount > 0
            SA.Enable = 'on';
        end
                
    case 1  % Object node selected
        set(handles.AddButton,'Enable','on','String','Add Object')
        set(handles.RemoveButton,'Enable','on','String','Remove Object')
                
    case 2  % Static/Dynamic node selected
        set(handles.AddButton,'Enable','on','String',addFileString())
        set(handles.RemoveButton,'Enable','off','String','Remove File')
                
    case 3  % File node selected
        set(handles.AddButton,'Enable','on','String',addFileString())
        set(handles.RemoveButton,'Enable','on','String','Remove File')
                
end %switch
set(handles.SmartAddButton,'Enable',SA.Enable)


    function aStr = addFileString()
        aStr = 'Add File(s)';
        node = first(tree.getSelectedNodes);        
        if node.getLevel == 3
            node = node.getParent;
        end
        cNodeStr = getCatName(node);    % 'Static' | 'Dynamic'
        if isequal( 'static', cNodeStr ) && ...
                node.getChildCount > 0
            aStr = 'Change File';
        end        
    end

end %btnStates()


% ------------------------------------------------------------------------
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
        
    case handles.OkButton
        
        models = getModels( handles );
        
        if ~isempty(models)
            % Check on dynamic files:
            %   - Dynamic list can be empty
            %   - But if it's non-empty, it must have the same number as every
            %       other model
            ndf = cellfun(@numel,{models.LoRes});  % full set
            ndf(ndf==0) = [];      % ignore empty lists
            
            if isempty(ndf) || all( ndf(1)==ndf )  % all must have the same number
                handles.output = models;
            else
                wstr = sprintf(...
                    ['The number of dynamic files is not consistent across ',...
                    'each of the models.  The set of dynamic files can be ',...
                    'empty for any given model, but if files are loaded ',...
                    'then all models must have the same number of files.',...
                    '\n\nCorrect the files and try again.\n']);
                warndlg(wstr,'Inconsistent number of dynamic files','modal')
                return
            end
        end
        
        % Ok, we're good to return the ouput:
        handles.output = models;
        
    otherwise
        return
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

end %figure1_CloseRequestFcn


% ------------------------------------------------------------------------
function AddButton_Cbk(hObject, eventdata)
%ADDBUTTON_CBK Callback for multi-funtional "Add" button for adding objects
%or files
%
handles = guidata(hObject);

tree = handles.tree;

% Node handles:
parent = first(tree.getSelectedNodes);  % Parent node    
rootNode = tree.getModel().getRoot();   % Root node

% Switch depending on what action is required:
switch lower(get(hObject,'String'))
    
    case 'add object'
        name = sprintf('Object %d',rootNode.getChildCount+1);
        name = inputdlg2('Name:','New Object',1,{name});
        
        % DRAWNOW seems to be required to prevent a hang on Matlab 2010b,
        % Windows XP 32:
        drawnow     
        
        if isempty(name) || isempty( strrep(name{1},' ', '') ) || ...
                nameInUse(tree,name)
            return
        else
            newObject(tree,rootNode,name)
        end
        
    case {'add file', 'add file(s)', 'change file'}
        % "Add File" can be enabled on a file or on a Static/Dynamic
        % container:
        if parent.getLevel == 3         % If a file is selected
            parent = parent.getParent;  % Get its parent
        end
        handles = addFiles(handles,tree,parent);
        if isempty(handles)
            return
        end
        
    otherwise
        error('unhandled button')
        
end

% Update handles:
guidata(hObject,handles)

% Focus tree:
handles.jtree.requestFocus


end %AddButton_Cbk()


% ------------------------------------------------------------------------
function SmartAddButton_Cbk(hObject, eventdata)
%SMARTADDBUTTON_CBK Callback for the Smart add files button.  This action
% allows multi-selection in a uigetfile dialog, then distributes those
% files to the appropriate Objects that exist.  If files do match object
% and category names, they are not dealt out and a list is displayed of the
% files that failed.



handles = guidata(hObject);
treeModel = handles.tree.getModel;
root = treeModel.getRoot();   % Root node

persistent fspec
if isempty( fspec )
    fspec = {...
        '*.mat', 'Mat files (*.mat)';
        '*.txt','Text files (*.txt)';
        '*.xyz','Text files (*.xyz)';
        '*.*',  'All Files (*.*)'};
end

[filename, pathname] = uigetfile( ...
    fspec, 'Select files to automatically associate', ...
    handles.userPath,...
    'MultiSelect', 'on');
if isequal(filename,0)
    return
elseif ~iscell(filename),
    filename = {filename};
else
    filename = filename(:);   % For convenience
end

% Most of the time in the following loop is taken up loading the data from
% the *.mat files using Cloud.Load()
% If we could pull that out here and objects into the loop directly, we
% could do the tree manipulation/node adding in one hit at the end and the
% user wouldn't see all the files loading individually.  But that's a fair
% bit more work...
objnames = lower( getChildNames(root) );
no = numel(objnames);
for oj = 1:no
    nf = numel(filename);
    inds = (1:nf)';
    id = ~cellfun(@isempty, strfind(filename,objnames{oj}));
    if any(id)
        objNode = root.getChildAt(oj-1);
        ojfiles = filename(id);
        static  = ~cellfun(@isempty,strfind(ojfiles,'static'));
        dynamic = ~cellfun(@isempty,strfind(ojfiles,'dynamic'));
        % Valid soution for static is to have only one file:
        if sum(static) > 1
            static(:) = false;  % Reset - don't add any of these files
        elseif sum(static) == 1
            % Add static files
            parent = getChildNodeByName(objNode,'static');
            handles = addStaticFileNode(handles,parent,[pathname ojfiles{static}]);
        end
        % Add dynamic files:
        if any(dynamic)
            parent = getChildNodeByName(objNode,'dynamic');
            handles = addDynamicFileNodes(handles,parent,pathname,ojfiles(dynamic));
        end
        done = static | dynamic;
        indsj = inds(id);
        filename(indsj(done)) = [];
    end
end

% Warn of files not assigned:
if ~isempty(filename)
    msg = ['The following files either did not match the objects in the list, or ',...
        'do not conform to the naming convention.  You will need to assign them ',...
        'manually:'];
    msg = [{msg}; ' '; filename(:)];
    msgbox(msg,'Files not assigned','help','modal');
end

% Update:
handles.userPath = pathname;
guidata(hObject,handles)

end %SmartAddButton_Cbk()

% ------------------------------------------------------------------------
function childNode = getChildNodeByName(parentNode,childName)
% Get the category node specified by 
names = getChildNames(parentNode);
ind = find(strcmpi(names,childName),1,'first');
childNode = parentNode.getChildAt(ind-1);
end %getChildNodeByName()

% ------------------------------------------------------------------------
function handles = addStaticFileNode(handles,parent,filepath)
% Add a file node to a "Static" cateogry node of object PARENT.

% Delete old file node, if one exists:
if parent.getChildCount > 0
    handles.tree.setSelectedNode( parent.getFirstChild )
    handles = removeSelectedNodes(handles);
end

% Then add a new file node:
handles = newFileNode(handles,handles.tree.getModel,parent,0,filepath);
end

% ------------------------------------------------------------------------
function handles = addDynamicFileNodes(handles,parent,pathname,filenames)
% Add files to the "Dynamic" cateogry node of object PARENT.  We insert the
% files listed in cell array FILENAMES into the list of objects that
% already exists (inserted based on filename)
%
% Now we insert these into the list that is already there, and at
% the same time add them to the handles.cList structure
assert(iscell(filenames),'FILENAMES must be a cell array of strings')
[filename,psns] = sortForInsert(parent,filenames);
for j = 1:numel(filename)
    handles = newFileNode(handles,handles.tree.getModel,parent,psns(j),[pathname,filenames{j}]);
end

    % ----------------------------
    function [flist,idx] = sortForInsert(parent,flist)
        % Inserts must be applied sequentially, and assumes the children of
        % PARENT are already sorted
        flist = sort_nat(flist);
        
        % Build cell array of existing node names:
        nk = parent.getChildCount;
        nodelist = cell(nk,1);
        for k = 1:nk
            nodelist{k} = char(parent.getChildAt(k-1).getName);
        end
        
        % Sequentially put each file into the list and get its index
        nf = numel(flist);
        idx = zeros(nf,1);
        for f = 1:nf
            [~,idxf] = sort_nat([nodelist;flist(f)]);
            idx(f) = find(idxf == numel(idxf)) + f - 2;
        end
    end %sortForInsert()

end
        
% ------------------------------------------------------------------------
function handles = addFiles(handles,tree,parent)
% Propmt the user to add files, then add them to the tree

treeModel = tree.getModel;

persistent fspec
if isempty( fspec )
    fspec = {...
        '*.mat', 'Mat files (*.mat)';
        '*.txt','Text files (*.txt)';
        '*.xyz','Text files (*.xyz)';
        '*.*',  'All Files (*.*)'};
end

objName = char(parent.getParent.getName);

% Change loading method based on Static/Dynamic category:
switch getCatName(parent) % 'Static' | 'Dynamic'
    
    case 'Static'
        % Only one 'static' file is permitted per object, so if a file
        % already exists, we first delete it, and then create a new one
        % from the user's selection.
        [filename, pathname, fidx] = uigetfile( ...
            fspec, [objName, ': Select high res static point cloud'], ...
            handles.userPath,...
            'MultiSelect', 'off');
        if isequal(filename,0)
            return
        end
        % Add the node:
        handles = addStaticFileNode(handles,parent,[pathname,filename]);
        
        
    case 'Dynamic'
        % The 'dynamic' list of files can be 0-N in size, and must remain
        % sorted by the file name.  We do this by inserting each new file
        % into the list in the correct location.
        [filenames, pathname, fidx] = uigetfile( ...
            fspec, [objName, ': Select low res dynamic point cloud(s)'], ...
            handles.userPath,...
            'MultiSelect', 'on');
        if isequal(filenames,0)
            return
        elseif ~iscell(filenames), 
            filenames = {filenames};
        else
            filenames = filenames(:);   % For convenience
        end
        handles = addDynamicFileNodes(handles,parent,pathname,filenames);
                      
end %switch

% Check for problems with loading:
if isempty(handles)
    return
end

updateCatName(tree,parent) % Set the number

% Remember the filter selection by re-sorting:
idx = 1:size(fspec,1);
idx(fidx) = [];
fspec = fspec( [fidx idx], :);

% Remember the path
handles.userPath = pathname;%(1:find(pathname == filesep,2,'last'));

% Expand tree to show added children
tree.setSelectedNode( parent.getFirstChild );

% Select parent
tree.setSelectedNode( parent );    

end %addFiles()


% ------------------------------------------------------------------------
function handles = newFileNode(handles,treeModel,parent,psn,fileOrCloud,q,x,tagAppend)
% Usage:
%   newFileNode(handles,treeModel,parent,psn,fname)        % Load dynamic model(s) from file
%   newFileNode(handles,treeModel,parent,psn,cloud)        % Load Dynamic model from cloud
%   newFileNode(handles,treeModel,parent,psn,fname,q,x)    % Static models - store q & x 
%   newFileNode(handles,treeModel,parent,psn,cloud,q,x)    % Static models - store q & x 


if nargout < 1
    error('Handles must be pushed back to caller')
end

ISSTATIC = strcmpi('Static', getCatName(parent) );

% Handle source:
if ~isa(fileOrCloud,'char')
    % Use cloud object directly:
    cloud = fileOrCloud;
    [~,fname,ext] = fileparts(cloud.Path);
    fname = [fname ext];
else % file path passed in:
    [pname,fname,ext] = fileparts(fileOrCloud);
    pname = [pname filesep];
    fname = [fname ext];
    [cloud,p] = Cloud.Load([pname fname],fname);   % Load cloud structure
    nc = numel(cloud);
    % For ROIs of Version 2 and up, we can save all phases of a dynamic
    % tracing in the same file, so here numel(cloud) = #_of_phases.  So we
    % recursively call this function to create new files nodes for each
    % cloud that was loaded:
    if ISSTATIC && nc > 1
        errordlg('This file contained a set of dynamic data, not static data.',...
            'Static data required','modal');
        handles = [];   % flag problem
        return
    end
    for c = 1:nc
        % Set the tag appendix only for dynamic clouds:
        if ISSTATIC
            if exist('q','var') && exist('x','var')
                % Loading static cloud which has motion:
                handles = newFileNode(handles,treeModel,parent,psn+(c-1),cloud(c),q,x);
            else
                % Loading static cloud without motion:
                handles = newFileNode(handles,treeModel,parent,psn+(c-1),cloud(c));
            end
        else
            if nc == 1
                % Loading one dynamic cloud from one file:
                handles = newFileNode(handles,treeModel,parent,psn+(c-1),cloud(c));
            else
                % Loading multiple dynamic clouds from one file:
                ta = sprintf('--phase-%02d',p(c));
                handles = newFileNode(handles,treeModel,parent,psn+(c-1),cloud(c),[],[],ta);
            end
        end            
        
    end
    return
end


% Store in hanldes.cList by uuid:
uuid = fieldUUID();                         % uuid
handles.cList.(uuid) = cloud;               % new field with uuid

% And create new motion vectors or populate with inputs:
if ISSTATIC
    if ~exist('q','var') && ~exist('x','var')
        q = [];
        x = [];
    end
    handles.qList.(uuid) = q;
    handles.xList.(uuid) = x;
end

% Create the name string:
name = fname;
if exist('tagAppend','var') && ~isempty(tagAppend)
    name = [name tagAppend];
end
    
% Create node & tag it with a uuid for this cloud:
node = uitreenode('v0', 'file', name, [], 0);

% Get appropriate icon:
switch ext
    case '.mat'
        ico = getIcon('matfile');
    otherwise 
        ico = getIcon('page');
end
node.setIcon(ico)           % Set icon
node.setUserObject(uuid);   % Store uuid in UserObject

treeModel.insertNodeInto(node,parent,psn);
end %newFileNode()


% ------------------------------------------------------------------------
function varargout = newObject(tree,parent,name)
% Usage:
%   newObject(tree,parent,name)
%   node = newObject(tree,parent,name)

treeModel = tree.getModel;

% Get the currently selected node:
selNode = tree.getSelectedNodes;

% First add the object node:
objNode = uitreenode('v0','object', name, [], 0);
treeModel.insertNodeInto(objNode,parent,parent.getChildCount()); % tack onto the end

% Create the "Static" and "Dynamic" nodes:
sNode = uitreenode('v0','selected', 'Static',  [], 0);
dNode = uitreenode('v0','selected', 'Dynamic', [], 0);

% Add them as children to the object node:
treeModel.insertNodeInto(sNode,objNode,objNode.getChildCount());
treeModel.insertNodeInto(dNode,objNode,objNode.getChildCount());

% Change the embedded icon:
jImg = getIcon('checked');
sNode.setIcon(jImg);
dNode.setIcon(jImg);

% Expand to show added children
tree.setSelectedNode( sNode );

% If root node was selected at the beginning, restore selction to that;
% otherwise, select the newly created node:
if ~isempty(selNode) && selNode(1).isRoot
    tree.setSelectedNode( selNode(1) );
else
    tree.setSelectedNode( objNode );
end

% Select Parent:
%tree.setSelectedNode( parent );

% Set optional output
if nargout > 0
    varargout{1} = objNode;
end

end %newObject()


% ------------------------------------------------------------------------
function nodelist = findnode(node,varargin)
% Example:
%   nodelist = findnode(rootnode,'getLevel',3) % find all nodes at level 3

if numel(varargin) < 2
    return
end

% Initialise node list:
nodelist = {};

% Drill down into the tree
props = varargin(1:2:end);
vals  = varargin(2:2:end);
for j = 1:numel(props)                  % Evaluate all properties (methods)
    nc = node.getChildCount();
    for cj = 1:nc
        kid = node.getChildAt(cj-1);
        if kid.(props{j}) == vals{j}    % Check on method
            nodelist{end+1,1} = kid;    % If true, add to list
        end
        if kid.getChildCount > 0        % If has kids, drill into kids
            nodelist = [nodelist; findnode(kid,varargin{:})];
        end
    end
end

end %findnode()


% ------------------------------------------------------------------------
function str = getCatName(catNode)
str = strtok( char(catNode.getName()) );
end %getCatName()


% ------------------------------------------------------------------------
function updateCatName(tree,catNode,baseName)
%UPDATECATNAME Set or update the category name
%
% Usage:
%   updateCatName(tree,catNode,baseName)     % For uninitialised nodes
%   updateCatName(tree,catNode)              % For initialised nodes

if nargin < 3
    baseName = getCatName(catNode);
    if isempty(baseName)
        keyboard %error()
    end
end

% Get how many files there used to be:
cstr = textscan(char(catNode.getName()),'%s (%d)');

nfiles = catNode.getChildCount;

if nfiles == 0
    %Name:  'Static' | 'Dynamic'
    name = sprintf('%s',baseName);  
else
    %Name:  'Static (1)' | 'Dynamic (5)'
    name = sprintf('%s (%d)',baseName,nfiles);
end
catNode.setName( name )

% Refresh node space like this:
refreshNodeSpace(catNode)

end %updateCatName()
 

% ------------------------------------------------------------------------
function refreshNodeSpace(node)
%REFRESHNODESPACE Refresh the horizontal space required for the node name
% The argument NODE is the node to refresh.
if isempty( node )
    return
end

handles = guidata( node(1).getRoot.getUserObject );
if ~all( isfield(handles,{'tree','jtree'}) )
    % Tree creation - fields not available yet
    return
end
tree  = handles.tree; 
jtree = handles.jtree;

% Refresh node space like this:
try % Regular method using nodeChanged:
    % This fails when we're building the tree for the first time, and I
    % haven't figure out why yet.
    tree.getModel.nodeChanged(node)
    
catch %#ok<CTCH>
    % So if that failed, we do a manual process:
    
    % Remember current selection & Multi-select state:
    drawnow
    selNode = first(tree.getSelectedNodes);
    ms = tree.isMultipleSelectionEnabled;
    
    % Set selection to node to refresh:
    tree.setMultipleSelectionEnabled( true )
    if numel(node) == 1
        tree.setSelectedNode( node )
    else
        tree.setSelectedNodes( node )
    end
    
    % Re-draw so jtree can pick up on the currently selected row:
    drawnow
    
    % Refreshing the horizontal space for the name
    % requires list collapse/expansion at this level:
    sr = jtree.getSelectionRows;
    for j = 1:numel(sr)
        if jtree.isCollapsed( sr(j) )
            jtree.expandRow( sr(j) );
            jtree.collapseRow( sr(j) );
        else
            jtree.collapseRow( sr(j) );
            jtree.expandRow( sr(j) );
        end
    end
    
    % Restore original selection & Multi-select state:
    tree.setSelectedNode( selNode );
    tree.setMultipleSelectionEnabled( ms );
end

end %refreshNodeSpace()


% ------------------------------------------------------------------------
function RemoveButton_Cbk(hObject, ~)
%REMOVEBUTTON_CBK Remove the currently selected item

handles = guidata(hObject);
handles = removeSelectedNodes(handles);

% Update handles:
guidata(hObject,handles)

% Set focus on tree instead of button:
handles.jtree.requestFocus

end %RemoveButton_Cbk()


% ------------------------------------------------------------------------
function handles = removeSelectedNodes(handles)
%REMOVESELECTEDNODE Remove the currently selected node from TREE and
% re-sets the selection.
tree = handles.tree;
nodes = tree.getSelectedNodes;

% Let's only allow removing of nodes that have the same parent:
for j = nodes.numel : -1 : 1
    parent(j) = nodes(j).getParent;
    same_parent(j) = parent(j) == parent(end);
end
if ~all(same_parent)
    errordlg(['When removing more than one object, the objects must ',...
        'all have the same parent. Change your selection and try again.'],...
        'Mixed selection not allowed','modal')
    return
end

% Root node cannot be removed (noting that if the root node is selected, it
% will be the only one of its kind, and by passing the test above it will
% be the only selected node):
if isRoot(first(nodes))
    return
end

% All the parents are the same, so collapse array to single:
parent = parent(1);

% Adjust selection:
%   Now that we have a node list, change the selection by selecting either
%   the next sibling (if there is one), or the previous sibling (if there
%   is one, or the parent (if there are no siblings):
nP = nodes(1).getPreviousSibling;
nN = nodes(end).getNextSibling;
if ~isempty( nN )
    nSel = nN;
elseif ~isempty( nP )
    nSel = nP;
else
    nSel = parent;
end
tree.setSelectedNode( nSel );

% Then remove the nodes:
for j = 1:numel(nodes)
    handles = removeThisNode(handles,nodes(j));
end

% Update the file count in category names:
if parent.getLevel == 2
    updateCatName(tree,parent)
end

end %removeSelectedNodes()

% ------------------------------------------------------------------------
function handles = removeThisNode(handles,node)
%REMOVETHISNODE Remove the specified node.
% No selection modification or other checking is done.

treeModel = handles.tree.getModel;

% If it's a 'file' node (at level 3), first remove the cloud
if node.getLevel == 3
    uuid = node.getUserObject();                    % Get cloud uuid
    handles.cList = rmfield(handles.cList, uuid);   % Remove the cloud
end

% Now we can delete the node
treeModel.removeNodeFromParent( node );         % Delete the node

end

% ------------------------------------------------------------------------
function ClearButton_Cbk(hObject, ~)

handles = guidata(hObject);
tree = handles.tree;
jtree = handles.jtree;

% Get root node:
rootnode = tree.getModel().getRoot();

% If file nodes exist, store them for deletion:
lv = 3;
nodes = findnode(rootnode,'getLevel',lv);

% Otherwise, store Object nodes for deletion
if isempty(nodes)
    lv = 1;
    nodes = findnode(rootnode,'getLevel',lv);
end
if isempty(nodes)
    return
end

%tree.setSelectedNode( rootnode );
for j = numel(nodes):-1:1
    handles = removeThisNode(handles,nodes{j});
end

% Refresh the category text
if lv == 3
    nodes = findnode(rootnode,'getLevel',2);
    for j = 1:numel(nodes)
        updateCatName(tree,nodes{j})
    end
end

guidata(hObject,handles)

end %ClearButton_Cbk()


% ------------------------------------------------------------------------
function OkButton_Cbk(hObject, ~)

% Close the figure:
figure1_CloseRequestFcn(hObject, [], guidata(hObject) )

end %OkButton_Cbk()


% ------------------------------------------------------------------------
function Models = getModels(handles)

tree = handles.tree;

% Get root node:
root = tree.getModel().getRoot();

% Pre-allocate space for Models:
% See also Session.m:
nm = root.getChildCount();
b = Bone();
Models = b(ones(1,nm)); % Faster than repmat

for mj = 1:nm
    mdlNode = root.getChildAt(mj-1); % 0-based indexing
    Models(mj).Tag = char(mdlNode.getName());
    
    for k = 1:mdlNode.getChildCount     % Step through the kinds / categories of files
        cNode = mdlNode.getChildAt(k-1);
        switch getCatName(cNode)    % 'Static' | 'Dynamic'
            case 'Static'
                [hr,q,x] = cloudStructFromNodeChildren( handles, cNode );
                Models(mj).HiRes = hr;
                Models(mj).q = q;
                Models(mj).x = x;
            case 'Dynamic'
                Models(mj).LoRes = cloudStructFromNodeChildren( handles, cNode );
        end
    end
end

% Now remove any empty models - ie, ones which have no static files, and no
% dynamic files: 
emT = cellfun(@isempty,{Models(:).LoRes}) & ...
    cellfun(@isempty,{Models(:).HiRes});
Models(emT) = [];

end %tree2models()


% ------------------------------------------------------------------------
function varargout = cloudStructFromNodeChildren( handles, parent )
% Return a cell list of file names from the child nodes of PARENT.  UUIDs
% are stored in the child node 'UserObject' property.

q = []; 
x = [];
cld = Cloud.empty();            % In case nf==0
nf = parent.getChildCount();            % Number of child nodes
for f = 1:nf                            % Cycle through child nodes    
    fnode = parent.getChildAt( f-1 );   % 0-based index
    uuid = fnode.getUserObject();       % UUID stored in UserObject
    cld(f,1) = handles.cList.(uuid);     % Get the cloud from the list
  
    % Now set the visiblity:
    switch parent.getValue
        case 'selected'
            vis = 1;
        case 'unselected'
            vis = 0;
    end
    cld(f,1).Visible = vis;
end

% Get (q,x) for Static files:
if nf == 1 && nargout == 3
    q = handles.qList.(uuid);
    x = handles.xList.(uuid);
end

if nargout >= 1
    varargout = {cld};
end
if nargout == 3
    varargout(2:3) = {q,x};
end

end %fileListFromNodeChildren()


% ------------------------------------------------------------------------
function CancelButton_Cbk(hObject, ~)

% Close the figure:
figure1_CloseRequestFcn(hObject, [], guidata(hObject))

end %CancelButton_Cbk()


% ------------------------------------------------------------------------
function a = first(vec)
%FIRST Select first element in array (usually array of nodes), or if empty
% return empty
% This is to circumvent the problem of doing the call vec(1) on an empty
% vector.
a = [];
if isempty(vec)
    return
end
if iscell(vec)
    a = vec{1};
else
    a = vec(1);
end
end %first()


% ------------------------------------------------------------------------
function uuid = fieldUUID()
% Create UUID which can be used as a field name.  This means:
%   - Cannot begin with a numeric
%   - Cannot include hyphens.
prefix = char(randi([97 122])); % random letter prefix
% prepend to avoid starting with a number:
uuid = [prefix char(java.util.UUID.randomUUID())];

% Remove hyphens:
uuid = strrep(uuid,'-','');

end %fieldUUID()


% ------------------------------------------------------------------------
function javaImage = getIcon(opt)

% Matlab icons can be found here:
%   [matlabroot '/toolbox/matlab/icons/']
switch opt
    case 'checked'
        I = uint8([...
            83,83,83,83,83,83,83,83,83,83,83,83,66,10,20
            83,83,83,83,83,83,83,83,83,83,83,82,11, 0, 9
            81,51,38,32,32,32,32,32,32,32,32,18, 0, 0,44
            49,41,60,78,79,80,80,80,80,79,60, 1, 0,14,83
            37,56,74,76,76,77,77,77,77,76,13, 0, 7,33,83
            29,69,69,69,70,71,71,71,71,34, 0, 0,39,29,83
            28,65,66,16, 3,40,67,67,54, 6, 0,17,69,28,83
            27,61,62, 8, 0, 5,45,64,15, 0, 4,57,61,27,83
            26,57,58,44, 3, 0, 9,31, 0, 0,35,60,57,26,83
            25,47,48,53,38, 0, 0, 0, 0,12,56,48,47,25,83
            24,50,51,52,57,30, 0, 0, 2,47,52,51,50,24,83
            22,55,56,57,57,65,23, 2,36,59,57,56,55,22,83
            21,59,61,62,63,63,71,75,68,63,62,61,59,21,83
            29,46,66,68,69,69,69,69,69,69,68,66,46,29,83
            43,34,51,72,73,74,74,74,74,73,72,51,34,43,83
            75,42,27,19,19,19,19,19,19,19,19,27,42,75,83]);
        map = [...
            0.09804, 0.09804, 0.09804
            0.11373, 0.11373, 0.11373
            0.11765, 0.11765, 0.11765
            0.12549, 0.12549, 0.12549
            0.13725, 0.13725, 0.13725
            0.14118, 0.14118, 0.14118
            0.14510, 0.14510, 0.14510
            0.15294, 0.15294, 0.15294
            0.20000, 0.20000, 0.20000
            0.21961, 0.21961, 0.21961
            0.24314, 0.24314, 0.24314
            0.26275, 0.26275, 0.26275
            0.27843, 0.27843, 0.27843
            0.29020, 0.29020, 0.29020
            0.31373, 0.31373, 0.31373
            0.32549, 0.32549, 0.32549
            0.32941, 0.32941, 0.32941
            0.33725, 0.33725, 0.33725
            0.36078, 0.36078, 0.36078
            0.38824, 0.38824, 0.38824
            0.41176, 0.41176, 0.41176
            0.41569, 0.41569, 0.41569
            0.42745, 0.42745, 0.42745
            0.43529, 0.43529, 0.43529
            0.44314, 0.44314, 0.44314
            0.45490, 0.45490, 0.45490
            0.47059, 0.47059, 0.47059
            0.48235, 0.48235, 0.48235
            0.49020, 0.49020, 0.49020
            0.49804, 0.49804, 0.49804
            0.50588, 0.50588, 0.50588
            0.52941, 0.52941, 0.52941
            0.53333, 0.53333, 0.53333
            0.54902, 0.54902, 0.54902
            0.56471, 0.56471, 0.56471
            0.56863, 0.56863, 0.56863
            0.58431, 0.58431, 0.58431
            0.58824, 0.58824, 0.58824
            0.60784, 0.60784, 0.60784
            0.62353, 0.62353, 0.62353
            0.63529, 0.63529, 0.63529
            0.65490, 0.65490, 0.65490
            0.66275, 0.66275, 0.66275
            0.67451, 0.67451, 0.67451
            0.72157, 0.72157, 0.72157
            0.72549, 0.72549, 0.72549
            0.73725, 0.73725, 0.73725
            0.74118, 0.74118, 0.74118
            0.74902, 0.74902, 0.74902
            0.75686, 0.75686, 0.75686
            0.76078, 0.76078, 0.76078
            0.76471, 0.76471, 0.76471
            0.76863, 0.76863, 0.76863
            0.77255, 0.77255, 0.77255
            0.78039, 0.78039, 0.78039
            0.78824, 0.78824, 0.78824
            0.79608, 0.79608, 0.79608
            0.80392, 0.80392, 0.80392
            0.81176, 0.81176, 0.81176
            0.81569, 0.81569, 0.81569
            0.81961, 0.81961, 0.81961
            0.82353, 0.82353, 0.82353
            0.83137, 0.83137, 0.83137
            0.83529, 0.83529, 0.83529
            0.83922, 0.83922, 0.83922
            0.84314, 0.84314, 0.84314
            0.84706, 0.84706, 0.84706
            0.85882, 0.85882, 0.85882
            0.86275, 0.86275, 0.86275
            0.87059, 0.87059, 0.87059
            0.87451, 0.87451, 0.87451
            0.87843, 0.87843, 0.87843
            0.88627, 0.88627, 0.88627
            0.89804, 0.89804, 0.89804
            0.90196, 0.90196, 0.90196
            0.90588, 0.90588, 0.90588
            0.90980, 0.90980, 0.90980
            0.91373, 0.91373, 0.91373
            0.93725, 0.93725, 0.93725
            0.94118, 0.94118, 0.94118
            0.94510, 0.94510, 0.94510
            0.94902, 0.94902, 0.94902
            0.98824, 0.98824, 0.98824
            1.00000, 1.00000, 1.00000];
        
    case 'unchecked'
        I = uint8([...
            28,28,28,28,28,28,28,28,28,28,28,28,28,28,28
            28,28,28,28,28,28,28,28,28,28,28,28,28,28,28
            24,11, 5, 2, 2, 2, 2, 2, 2, 2, 2, 5,11,24,28
            10, 7,14,28,28,28,28,28,28,28,28,14, 7,10,28
             4,14,28,28,28,28,28,28,28,28,28,28,14, 4,28
             1,28,28,28,28,28,28,28,28,28,28,28,28, 1,28
             0,27,27,27,27,27,27,27,27,27,27,27,27, 0,28
             0,26,25,25,25,25,25,25,25,25,25,25,26, 0,28
             0,25,24,23,23,23,23,23,23,23,23,24,25, 0,28
             0,21,16,15,15,15,15,15,15,15,15,16,21, 0,28
             0,21,16,15,15,15,15,15,15,15,15,16,21, 0,28
             0,22,17,15,15,15,15,15,15,15,15,17,22, 0,28
             0,23,18,16,16,16,16,16,16,16,16,18,23, 0,28
             3,13,19,18,17,17,17,17,17,17,18,19,13, 3,28
             9, 6,12,21,20,20,20,20,20,20,21,12, 6, 9,28
            17, 8, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 8,17,28]);
        map = [...
            0.61176, 0.61176, 0.61176
            0.61961, 0.61961, 0.61961
            0.65098, 0.65098, 0.65098
            0.66667, 0.66667, 0.66667
            0.68627, 0.68627, 0.68627
            0.70980, 0.70980, 0.70980
            0.72941, 0.72941, 0.72941
            0.74902, 0.74902, 0.74902
            0.77647, 0.77647, 0.77647
            0.78431, 0.78431, 0.78431
            0.81961, 0.81961, 0.81961
            0.82353, 0.82353, 0.82353
            0.87059, 0.87059, 0.87059
            0.87451, 0.87451, 0.87451
            0.90588, 0.90588, 0.90588
            0.92549, 0.92549, 0.92549
            0.92941, 0.92941, 0.92941
            0.93333, 0.93333, 0.93333
            0.93725, 0.93725, 0.93725
            0.94510, 0.94510, 0.94510
            0.94902, 0.94902, 0.94902
            0.95294, 0.95294, 0.95294
            0.95686, 0.95686, 0.95686
            0.96078, 0.96078, 0.96078
            0.96471, 0.96471, 0.96471
            0.97647, 0.97647, 0.97647
            0.98431, 0.98431, 0.98431
            0.99216, 0.99216, 0.99216
            1.00000, 1.00000, 1.00000
            1.00000, 1.00000, 1.00000
            1.00000, 1.00000, 1.00000
            1.00000, 1.00000, 1.00000];
        
    case 'page'
        I = uint8([...
            47,44,41,41,40,31,20,30,36,42,59,81,81,81,81,81
            34,52,54,54,50,44,40,22,39,57,28,53,81,81,81,81
            37,54,56,55,55,52,45,39,21,79,57,25,52,81,81,81
            37,57,58,57,57,56,52,44, 6,79,79,57,27,62,81,81
            37,60,61,61,60,59,58,51,17, 6, 4,10,18,46,81,81
            37,63,63,46,54,65,48,60,53,42,30,29,11,32,81,81
            37,66,67,19,18,58,16,13,25,28,48,41,29,27,81,81
            37,68,69,67,47,64,63,67,62,68,71,49,41,23,81,81
            37,70,70,28, 6,24, 8,11, 3,26,37,59,48,18,81,81
            37,70,72,75,56,68,68,54,76,64,55,66,60,18,81,81
            37,71,74,43,15, 9, 5,14,38,12, 5,69,66,18,81,81
            37,76,75,67,47,64,63,67,62,68,71,71,69,18,81,81
            37,76,78,28, 7,25, 8,11, 3,26,37,74,71,18,81,81
            37,78,80,75,56,68,68,54,76,64,55,76,73,18,81,81
            35,80,79,79,79,79,78,78,77,77,76,76,75,16,81,81
            33, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2,28,81,81]);
        map = [...
            0.11765, 0.11765, 0.11765
            0.14510, 0.14510, 0.14510
            0.25098, 0.25098, 0.25098
            0.54902, 0.54902, 0.54902
            0.56471, 0.56471, 0.56471
            0.58431, 0.58431, 0.58431
            0.60392, 0.60392, 0.60392
            0.60784, 0.60784, 0.60784
            0.61961, 0.61961, 0.61961
            0.63922, 0.63922, 0.63922
            0.64706, 0.64706, 0.64706
            0.65098, 0.65098, 0.65098
            0.65490, 0.65490, 0.65490
            0.66275, 0.66275, 0.66275
            0.67059, 0.67059, 0.67059
            0.67451, 0.67451, 0.67451
            0.68627, 0.68627, 0.68627
            0.69020, 0.69020, 0.69020
            0.69412, 0.69412, 0.69412
            0.69804, 0.69804, 0.69804
            0.70196, 0.70196, 0.70196
            0.70588, 0.70588, 0.70588
            0.70980, 0.70980, 0.70980
            0.71373, 0.71373, 0.71373
            0.71765, 0.71765, 0.71765
            0.72157, 0.72157, 0.72157
            0.72549, 0.72549, 0.72549
            0.73333, 0.73333, 0.73333
            0.73725, 0.73725, 0.73725
            0.74118, 0.74118, 0.74118
            0.74510, 0.74510, 0.74510
            0.74902, 0.74902, 0.74902
            0.75686, 0.75686, 0.75686
            0.76078, 0.76078, 0.76078
            0.76471, 0.76471, 0.76471
            0.76863, 0.76863, 0.76863
            0.77255, 0.77255, 0.77255
            0.77647, 0.77647, 0.77647
            0.78431, 0.78431, 0.78431
            0.79608, 0.79608, 0.79608
            0.80000, 0.80000, 0.80000
            0.81569, 0.81569, 0.81569
            0.81961, 0.81961, 0.81961
            0.83529, 0.83529, 0.83529
            0.83922, 0.83922, 0.83922
            0.84314, 0.84314, 0.84314
            0.84706, 0.84706, 0.84706
            0.85490, 0.85490, 0.85490
            0.85882, 0.85882, 0.85882
            0.86667, 0.86667, 0.86667
            0.87059, 0.87059, 0.87059
            0.87451, 0.87451, 0.87451
            0.87843, 0.87843, 0.87843
            0.88235, 0.88235, 0.88235
            0.88627, 0.88627, 0.88627
            0.89020, 0.89020, 0.89020
            0.89412, 0.89412, 0.89412
            0.89804, 0.89804, 0.89804
            0.90196, 0.90196, 0.90196
            0.90588, 0.90588, 0.90588
            0.90980, 0.90980, 0.90980
            0.91373, 0.91373, 0.91373
            0.92157, 0.92157, 0.92157
            0.92549, 0.92549, 0.92549
            0.92941, 0.92941, 0.92941
            0.93333, 0.93333, 0.93333
            0.93725, 0.93725, 0.93725
            0.94118, 0.94118, 0.94118
            0.94510, 0.94510, 0.94510
            0.94902, 0.94902, 0.94902
            0.95294, 0.95294, 0.95294
            0.95686, 0.95686, 0.95686
            0.96078, 0.96078, 0.96078
            0.96471, 0.96078, 0.96078
            0.96471, 0.96471, 0.96471
            0.97255, 0.97255, 0.97255
            0.97647, 0.97647, 0.97647
            0.98039, 0.98039, 0.98039
            0.98431, 0.98431, 0.98431
            0.98824, 0.98824, 0.98824
            0.99216, 0.99216, 0.99216
            1.00000, 1.00000, 1.00000];
    case 'matfile'
        I = uint8([...
            74,74,33,24,24,24,24,24,24,24,24,24,74,74,74,74
            74,74,33,69,13,14,45,72,71,69,67,24,24,74,74,74
            74,74,33,11, 3,19,18,70,69,69,69,24,53,24,74,74
            65,51,10, 0, 2,16, 7,59,69,69,69,24,24,24,24,74
            29,12, 9, 1, 5,15, 8,31,64,62,61,60,50,48,24,56
            74,49,17, 4, 6,23,47,20,27,27,27,27,27,57,24,56
            74,74,33,21,22,54,46,46,36,37,39,43,27,57,24,56
            74,74,33,52,55,74,58,46,74,25,73,63,27,57,24,56
            74,74,33,57,27,28,34,37,41,38,40,44,27,57,24,56
            74,74,33,57,27,74,74,37,74,38,74,68,27,57,24,56
            74,74,33,57,27,28,28,26,38,30,32,35,27,57,24,56
            74,74,33,57,27,74,74,42,74,43,74,66,27,57,24,56
            74,74,33,57,27,27,27,27,27,27,27,27,27,57,24,56
            74,74,33,57,57,57,57,57,57,57,57,57,57,57,24,56
            74,74,33,24,24,24,24,24,24,24,24,24,24,24,24,56
            74,74,74,56,56,56,56,56,56,56,56,56,56,56,56,56]);
        map = [...
            0.33333, 0.20392, 0.29412
            0.45098, 0.19608, 0.23529
            0.60392, 0.19608, 0.14510
            0.42745, 0.25098, 0.30980
            0.66667, 0.23922, 0.16078
            0.72157, 0.28235, 0.18824
            0.85882, 0.42353, 0.24706
            0.87843, 0.41569, 0.27843
            0.90588, 0.45882, 0.26275
            0.35686, 0.58039, 0.76078
            0.45490, 0.57255, 0.70980
            0.46667, 0.61569, 0.72941
            0.33333, 0.66667, 0.86275
            0.56471, 0.61961, 0.70196
            0.86275, 0.58431, 0.25490
            0.94118, 0.56471, 0.28627
            0.93333, 0.57255, 0.28235
            0.65882, 0.62353, 0.66667
            0.89804, 0.56863, 0.50980
            0.91765, 0.60784, 0.26667
            0.85098, 0.62353, 0.60784
            0.92941, 0.65098, 0.38039
            0.93333, 0.65098, 0.41961
            0.92157, 0.65882, 0.47059
            0.65098, 0.72157, 0.79608
            0.75294, 0.75294, 0.76471
            0.75294, 0.75686, 0.76471
            0.72549, 0.76078, 0.81961
            0.75686, 0.75686, 0.76863
            0.65098, 0.78039, 0.85490
            0.76078, 0.76078, 0.76863
            0.90980, 0.72549, 0.70980
            0.76863, 0.76863, 0.77255
            0.72549, 0.77647, 0.83529
            0.77255, 0.77647, 0.79216
            0.78824, 0.77647, 0.77255
            0.78039, 0.78431, 0.79216
            0.78431, 0.78431, 0.79216
            0.78431, 0.78824, 0.79608
            0.78824, 0.79216, 0.79608
            0.79216, 0.79216, 0.80000
            0.79608, 0.80000, 0.81176
            0.80392, 0.80000, 0.79608
            0.80784, 0.80000, 0.79608
            0.81176, 0.80392, 0.79608
            0.93725, 0.78039, 0.76471
            0.81569, 0.81961, 0.82745
            0.87843, 0.81569, 0.78039
            0.79608, 0.83922, 0.87843
            0.80000, 0.84314, 0.87059
            0.80784, 0.84314, 0.88235
            0.75294, 0.87843, 0.94510
            0.89020, 0.85490, 0.77647
            0.83922, 0.87059, 0.90196
            0.88627, 0.87451, 0.86275
            0.89020, 0.88235, 0.86275
            0.90196, 0.92549, 0.94510
            0.94118, 0.95294, 0.96078
            0.94510, 0.95294, 0.96078
            0.97647, 0.94902, 0.94902
            0.95686, 0.96471, 0.97255
            0.96078, 0.96863, 0.97647
            0.96863, 0.97255, 0.97647
            0.96863, 0.97647, 0.97647
            0.97255, 0.97647, 0.98039
            0.97255, 0.98039, 0.98431
            0.98039, 0.98824, 0.98431
            0.98824, 0.98824, 0.99216
            0.98039, 0.99216, 0.99608
            0.99216, 0.99216, 0.99216
            0.99216, 0.99216, 0.99608
            0.99608, 0.99216, 0.99216
            0.99608, 0.99608, 0.99608
            0.99216, 1.00000, 1.00000
            1.00000, 1.00000, 1.00000];

    otherwise
        error('bad option')
end

javaImage = im2java(I,map);

end %getIcon()

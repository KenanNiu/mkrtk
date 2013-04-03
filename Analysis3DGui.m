function varargout = Analysis3DGui(varargin)
% ANALYSIS3DGUI MATLAB code for Analysis3DGui.fig
%      ANALYSIS3DGUI, by itself, creates a new ANALYSIS3DGUI or raises the existing
%      singleton*.
%
%      H = ANALYSIS3DGUI returns the handle to a new ANALYSIS3DGUI or the handle to
%      the existing singleton*.
%
%      ANALYSIS3DGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYSIS3DGUI.M with the given input arguments.
%
%      ANALYSIS3DGUI('Property','Value',...) creates a new ANALYSIS3DGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Analysis3DGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Analysis3DGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Analysis3DGui

% Last Modified by GUIDE v2.5 22-Oct-2012 15:10:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Analysis3DGui_OpeningFcn, ...
                   'gui_OutputFcn',  @Analysis3DGui_OutputFcn, ...
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


% --- Executes just before Analysis3DGui is made visible.
function Analysis3DGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Analysis3DGui (see VARARGIN)

% Choose default command line output for Analysis3DGui
handles.output = hObject;

% Bail out (but show figure) if inputs are incorrect:
if numel(varargin) == 0
    delete(handles.figure1)
    error([upper(mfilename) ' must be called with inputs.']);
end


% Here we convert the input structure to the form required by the
% listboxes.  This input form is the same as the output form, so see
% figure1_CloseRequestFcn for the output conversions and exporting the
% structure.

% Now, VARARGIN should contain the following:
%   #1  1-by-N cell array of object tags which are candidates for using in helical axis definitions 
%   #2  1-by-N cell array of object tags which are candidates for using in line of action definitions 
%   #3  1-by-N struct array of HelicalAxis definitions
%   #4  1-by-N struct array of LineOfAction definitions
%   #5  1-by-N strucy array of MomentArm definitions
%
% As follows:
%   varargin = { HaxCandTags, LoaCandTags, HelicalAxes, LinesOfAction, MomentArms }
    

% Manage input #1: Cell array of object tags which are candidates for being
%                  used in helical axis calculations:
handles.AxisCandidates = varargin{1}(:);    % Bone Names in column format
set(handles.HaxObj1Popup,'String',handles.AxisCandidates)   %\ Populate listboxes
set(handles.HaxObj2Popup,'String',handles.AxisCandidates)   %/ 

% Manage input #2: Cell array of object tags which are candidates for being
%                  used in line of action calculations:
handles.LineCandidates = varargin{2}(:);
set(handles.LoaObjPopup,'String',handles.LineCandidates)      

% Deal other inputs:
[hax,loa,marm] = varargin{3:end};

% Store for checking outputs:
handles.inputs.hax = hax;
handles.inputs.loa = loa;
handles.inputs.marm = marm;

% Write info into gui uicontrols:
populateList('Hax',hax,handles.HaxDefnListbox,handles.DisplayHaxCheckbox);
populateList('Loa',loa,handles.LoaDefnListbox,handles.DisplayLoaCheckbox);
populateList('Marm',marm,handles.MarmDefnListbox,handles.DisplayMarmCheckbox);

% Position gui over the top of calling gui
positionOver(handles.figure1,gcbf);

% Configure some button pictures:
load('icons.mat');
set(handles.NewHaxButton,'CData',icons.rightarrow)
set(handles.DelHaxButton,'CData',icons.delete_grey)
set(handles.NewLoaButton,'CData',icons.rightarrow)
set(handles.DelLoaButton,'CData',icons.delete_grey)
set(handles.NewMarmButton,'CData',icons.rightarrow)
set(handles.DelMarmButton,'CData',icons.delete_grey)

% Add some key definitions:
handles.dklist = {'delete','del','backspace'};  % Delete key options

% Add some tooltips
set(handles.HaxDefnListbox,'Tooltip',ttgen('HAx'))
set(handles.LoaDefnListbox,'Tooltip',ttgen('LoA'))
set(handles.MarmDefnListbox,'Tooltip',ttgen('Marm'))

% Listeners to control list selection (value):
addlistener(handles.HaxDefnListbox, 'String','PostSet',@selectionInBounds);
addlistener(handles.LoaDefnListbox, 'String','PostSet',@selectionInBounds);
addlistener(handles.MarmDefnListbox,'String','PostSet',@selectionInBounds);
addlistener(handles.MarmObj1Popup,  'String','PostSet',@selectionInBounds);
addlistener(handles.MarmObj2Popup,  'String','PostSet',@selectionInBounds);

% Listeners to keep moment arm stuff up to date:
addlistener(handles.HaxDefnListbox,'String','PostSet',@updateMarmPopups);
addlistener(handles.LoaDefnListbox,'String','PostSet',@updateMarmPopups);
addlistener(handles.HaxDefnListbox,'Value','PostSet',@updateMarmPopups);
addlistener(handles.LoaDefnListbox,'Value','PostSet',@updateMarmPopups);
updateMarmPopups(handles.HaxDefnListbox);
updateMarmPopups(handles.LoaDefnListbox);

% Listeners to keep the moment arm 
addlistener(handles.HaxDefnListbox,'Value','PostSet',@updateHaxPopups);
updateHaxPopups(handles.HaxDefnListbox);

% Focus "OK" button
uicontrol(handles.OkButton)

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = Analysis3DGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in HaxObj1Popup.
function HaxObj1Popup_Callback(hObject, eventdata, handles)
% hObject    handle to HaxObj1Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns HaxObj1Popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from HaxObj1Popup


% --- Executes during object creation, after setting all properties.
function HaxObj1Popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HaxObj1Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in HaxDefnListbox.
function HaxDefnListbox_Callback(hObject, eventdata, handles)
% hObject    handle to HaxDefnListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi( 'open', get(handles.figure1,'SelectionType'))
    userEditListItem(hObject)
end

% --- Executes during object creation, after setting all properties.
function HaxDefnListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HaxDefnListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on key press with focus on HaxDefnListbox and none of its controls.
function HaxDefnListbox_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to HaxDefnListbox (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if isempty(eventdata.Modifier) && ...
        any(strcmpi(eventdata.Key,handles.dklist))
        % Delete currently selected line item:
        runCallback(handles.DelHaxButton)
end


% --- Executes on selection change in HaxObj2Popup.
function HaxObj2Popup_Callback(hObject, eventdata, handles)
% hObject    handle to HaxObj2Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns HaxObj2Popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from HaxObj2Popup


% --- Executes during object creation, after setting all properties.
function HaxObj2Popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HaxObj2Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NewHaxButton.
function NewHaxButton_Callback(hObject, eventdata, handles)
% hObject    handle to NewHaxButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cstr = get(handles.HaxDefnListbox,'String');

if numel(cstr) >= 9
    warndlg('Only up to 9 axis definitions permitted','Max definitions permitted','modal')
    return
end

obj1 = current_item(handles.HaxObj1Popup);
obj2 = current_item(handles.HaxObj2Popup);

cstr = editList('add','Hax',cstr,obj1,obj2);

% Update list box:
set(handles.HaxDefnListbox,'String',cstr);



% --- Executes on button press in DelHaxButton.
function DelHaxButton_Callback(hObject, eventdata, handles)
% hObject    handle to DelHaxButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get listbox info
val = get(handles.HaxDefnListbox,'Value');
cstr = get(handles.HaxDefnListbox,'String');

% Edit the list and record the old and new information
old = html2cell(cstr);                      % Contents as a cell removing item
cstr = editList('remove','Hax',cstr,val);   % Remove item
new = html2cell(cstr);                      % Contents as a cell after item

% Now update the Moment Arm listbox using the changes:
newList = [new{1}(1:val-1); {[]}; new{1}(val:end)]; % Create blank space where item removed
oldList = old{1};                                   % New variable for clarity
MarmListUpdater(handles,oldList,newList)            % Update the Moment Arm list

% Update
set(handles.HaxDefnListbox,'String',cstr);

% Focus the listbox
uicontrol(handles.HaxDefnListbox)


% --- Executes on button press in DisplayHaxCheckbox.
function DisplayHaxCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayHaxCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DisplayHaxCheckbox


% --- Executes on selection change in LoaDefnListbox.
function LoaDefnListbox_Callback(hObject, eventdata, handles)
% hObject    handle to LoaDefnListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi( 'open', get(handles.figure1,'SelectionType'))
    userEditListItem(hObject)
end

% --- Executes during object creation, after setting all properties.
function LoaDefnListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoaDefnListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in LoaObjPopup.
function LoaObjPopup_Callback(hObject, eventdata, handles)
% hObject    handle to LoaObjPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LoaObjPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LoaObjPopup


% --- Executes during object creation, after setting all properties.
function LoaObjPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoaObjPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on key press with focus on LoaDefnListbox and none of its controls.
function LoaDefnListbox_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to LoaDefnListbox (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if isempty(eventdata.Modifier) && ...
        any(strcmpi(eventdata.Key,handles.dklist))
    % Delete currently selected line item:
    runCallback(handles.DelLoaButton);
end


% --- Executes on button press in NewLoaButton.
function NewLoaButton_Callback(hObject, eventdata, handles)
% hObject    handle to NewLoaButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cstr = get(handles.LoaDefnListbox,'String');

if numel(cstr) >= 9
    warndlg('Only up to 9 loading axes permitted','Max definitions','modal')
    return
end

obj = current_item(handles.LoaObjPopup);

cstr = editList('add','LoA',cstr,obj);

% Update listbox:
set(handles.LoaDefnListbox,'String',cstr);


% --- Executes on button press in DelLoaButton.
function DelLoaButton_Callback(hObject, eventdata, handles)
% hObject    handle to DelLoaButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get listbox info
val = get(handles.LoaDefnListbox,'Value');
cstr = get(handles.LoaDefnListbox,'String');

% Edit the list and record the old and new information
old = html2cell(cstr);                      % Contents as a cell removing item
cstr = editList('remove','loa',cstr,val);   % Remove item
new = html2cell(cstr);                      % Contents as a cell after item

% Now update the Moment Arm listbox using the changes:
newList = [new{1}(1:val-1); {[]}; new{1}(val:end)]; % Create blank space where item removed
oldList = old{1};                                   % New variable for clarity
MarmListUpdater(handles,oldList,newList)            % Update the Moment Arm list

% Update
set(handles.LoaDefnListbox,'String',cstr)

% Focus listbox:
uicontrol(handles.LoaDefnListbox)


% --- Executes on button press in DisplayLoaCheckbox.
function DisplayLoaCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayLoaCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DisplayLoaCheckbox


% --- Executes on selection change in MarmDefnListbox.
function MarmDefnListbox_Callback(hObject, eventdata, handles)
% hObject    handle to MarmDefnListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi( 'open', get(handles.figure1,'SelectionType'))
    userEditListItem(hObject)
end



% --- Executes during object creation, after setting all properties.
function MarmDefnListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MarmDefnListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on MarmDefnListbox and none of its controls.
function MarmDefnListbox_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to MarmDefnListbox (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if isempty(eventdata.Modifier) && ...
        any(strcmpi(eventdata.Key,handles.dklist))
    % Delete currently selected line item:
    runCallback(handles.DelMarmButton);
end


% --- Executes on button press in NewMarmButton.
function NewMarmButton_Callback(hObject, eventdata, handles)
% hObject    handle to NewMarmButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

obj1 = current_item(handles.HaxDefnListbox);
obj2 = current_item(handles.LoaDefnListbox);

% Check that the user has defined other geometry.
if isempty(obj1) || isempty(obj2)
    return
end

% Get string components:
parts1 = html2cell(obj1);
parts2 = html2cell(obj2);

% Add new component to list:
cstr = get(handles.MarmDefnListbox,'String');
cstr = editList('add','marm',cstr,parts1{1},parts2{1});

% Update list
set(handles.MarmDefnListbox,'String',cstr);


% --- Executes on button press in DelMarmButton.
function DelMarmButton_Callback(hObject, eventdata, handles)
% hObject    handle to DelMarmButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(handles.MarmDefnListbox,'Value');
cstr = get(handles.MarmDefnListbox,'String');
cstr = editList('remove','Marm',cstr,val);

% Update
set(handles.MarmDefnListbox,'String',cstr);


% --- Executes on button press in DisplayMarmCheckbox.
function DisplayMarmCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayMarmCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DisplayMarmCheckbox



% --- Executes on button press in ClearAllButton.
function ClearAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clear inputs:
f = fields(handles.inputs);
for j = 1:numel(f)
    handles.inputs.(f{j}) = [];
end
% Clear listboxes:
C ={};
set(handles.HaxDefnListbox,'String',C);
set(handles.LoaDefnListbox,'String',C);
set(handles.MarmDefnListbox,'String',C);
% Update data:
guidata(hObject,handles)


% --- Executes on button press in OkButton.
function OkButton_Callback(hObject, eventdata, handles)
% hObject    handle to OkButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure1_CloseRequestFcn(hObject, [], handles)


% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure1_CloseRequestFcn(hObject, [], handles)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Manage specific case: Window closed with title "X" button:
if isempty(handles)
    delete(hObject)
    return
end

% Switch the output depending on how the user closes the dialog:
switch hObject
    
    case handles.figure1
        % In this case, user closed using the window X button
        % OUTPUT: GUI figure handle (default)
        
    case handles.CancelButton
        % User pressed cancel.
        % OUTPUT: GUI figure handle (default)
        
    case handles.OkButton
        % User pressed ok:
        % OUTPUT: Geometry defitions from listboxes, with calculations
        % preserved from inputs
        %
        % **A note on outputs:**
        %
        % The neat thing about the way this output works is that if the
        % user deletes a definition, and then defines a new one using the
        % same objects, the data from the old one will be preserved on
        % output.  Not helpful if there's a problem with that original one,
        % but that shouldn't happen often.  If there is, the work-around
        % would be to delete it, close this GUI, then reopen and define the
        % new one.  Otherwise, preserving calculations is a good thing,
        % because they can be computationally expensive.
        
        try
        % Build output:
        hax = html2cell(get(handles.HaxDefnListbox,'String'));
        loa = html2cell(get(handles.LoaDefnListbox,'String'));
        marm = html2cell(get(handles.MarmDefnListbox,'String'));
        
        vis = [get(handles.DisplayHaxCheckbox,'Value'),...
            get(handles.DisplayLoaCheckbox,'Value'),...
            get(handles.DisplayMarmCheckbox,'Value')];
        
        % Here we build the ouput structure.  This output form will be the
        % same as the input form, so see the OpeningFcn as well for reading
        % the structure in.

        Helical = struct('Tag',{},...
            'Item1',{},...
            'Item2',{},...
            'Visible',{},...
            'Axis',{},...
            'Angle',{},...
            'Point',{},...
            'Slide',{});
        LineAction = struct('Tag',{},...
            'Item',{},...
            'Visible',{},...
            'Point',{},...
            'Vector',{});
        MomentArm = struct('Tag',{},...
            'Item1',{},...
            'Item2',{},...
            'Visible',{},...
            'Point1',{},...
            'Point2',{},...
            'Length',{});
        
        % Expand to size:
        [Helical(1:numel(hax{1})).Tag] = deal('');
        [LineAction(1:numel(loa{1})).Tag] = deal('');
        [MomentArm(1:numel(marm{1})).Tag] = deal('');
        
        % Helper function:
        getrow = @(defn,j)cellfun(@(x)x{j},defn,'UniformOutput',false);
        
        % Create / update Helical Axes:
        for j = 1:numel(hax{1})
            Helical(j) = configureOutputStruct(Helical(j),handles.inputs.hax,getrow(hax,j),vis(1));
        end
        
        % Create / update Lines of Action:
        for j = 1:numel(loa{1})
            LineAction(j) = configureOutputStruct(LineAction(j),handles.inputs.loa,getrow(loa,j),vis(2));
        end
        
        % Create / update Moment Arms:
        for j = 1:numel(marm{1})
            MomentArm(j) = configureOutputStruct(MomentArm(j),handles.inputs.marm,getrow(marm,j),vis(3));
        end
        
        catch ME
            % Around 22-Oct-12 we renamed some fields like 'Pt' -> 'Point'
            % and 'Vec' => Vector, so let's supply an informative dialog if
            % that looks to be the cause of the error:
            if isequal(ME.identifier,'MATLAB:heterogeneousStrucAssignment')
                errordlg(['Looks like you might have defined some of these '...
                    'with an earlier version.  Try clearing and re-defining them.'],...
                    'Old definitions?','modal')
            end
            rethrow(ME)
        end
        
        % Populate output field:
        handles.output = {Helical,LineAction,MomentArm};
        
        % Update gui
        guidata(hObject,handles)
        
end

% First, hide the figure so the user can't mess anything up:
set(handles.figure1,'Visible','off')

% Next, trigger the listener from invoking program, MRIMagic:
set(handles.figure1,'HitTest','off')

% Then when the listener completes, it will return here and delete the
% figure:
delete(handles.figure1)


% ------------------------------------------------------------------------
function new = configureOutputStruct(new,old,usrdefn,vis)
% Update the structure NEW from either an equivalent OLD structure or from
% a USRDEFN
% NEW and OLD are structures.
% USRDEFN is a cell like:
%   {'Axis_1', 'Tibia', 'Talus'}
%   {'Line_1', 'Achilles'}
if isempty(old)
    idx = [];
else
    if isfield(old,'Item')
        idx = strcmpi(usrdefn{2},{old.Item});
    elseif isfield(old,'Item1')
        idx = strcmpi(usrdefn{2},{old.Item1}) & strcmpi(usrdefn{3},{old.Item2});
    else
        error('Incorrect usage')
    end
end
if any(idx)
    % Found a matching input - copy it across to output:
    new = old(idx);
else
    % No matching inputs - define new:
    if isfield(new,'Item')
        new.Item = usrdefn{2};
    elseif isfield(new,'Item1')
        new.Item1 = usrdefn{2};
        new.Item2 = usrdefn{3};
    end
end
% In either case, update the Tag & Visibility
new.Tag = usrdefn{1};
new.Visible = vis(1);

            
% ------------------------------------------------------------------------
function populateList(opt,hstruct,listbox,checkbox)
% Populate the listbox specified by OPT with the data from HSTRUCT
%
% OPT is one of the options as per EDITLIST.
%
% HSTRUCT is the input/output structure of any of the axis definitions

if ~isempty(hstruct)
    
    % Build the string
    cstr = {};
    for j = 1:numel(hstruct)
        if isfield(hstruct,'Item2')
            cstr = editList('add',opt,cstr,hstruct(j).Item1,hstruct(j).Item2);
        else
            cstr = editList('add',opt,cstr,hstruct(j).Item);
        end
    end
    
    % Set the string:
    set(listbox,'String',cstr);
    
    % Set the visibility:
    set(checkbox,'Value',hstruct(1).Visible)
    
end


% ------------------------------------------------------------------------
function userEditListItem(hobj)
% Edit the currently selected item from the list
handles = guidata(hobj);
sel = get(hobj,'Value');
S = get(hobj,'String');
C = html2cell(S);

c = {C{1}{sel} C{2}{sel} C{3}{sel}};
n = sum(~cellfun(@isempty,c))-1;        % Number of popups

% Build interface:
bkgclr = get(gcbf,'Color');
hf = figure('Name','Edit property',...
    'Visible','off',...
    'Toolbar','none',...
    'Menubar','none',...
    'Color',bkgclr);

h = 20;     % uicontrol height
y = 20;     % base y-value
tw = 60;    % text width
bw = 50;    % button width
pw = 130;   % popup width
sep = 10;   % separation
x = 5;      % initial x value
w = tw;     % initial width
ht = uicontrol(hf,'Style','text',...
    'Position',[x y w h],...
    'HorizontalAlignment','right',...
    'BackgroundColor',bkgclr,...
    'FontWeight','bold',...
    'String',[c{1} ':']);
x = x+w+3;
w = pw;

% Create popup contents:
switch get(hobj,'Tag')
    case 'HaxDefnListbox'
        tag = 'Hax';
        str{1} = handles.AxisCandidates;
        str{2} = str{1};
    case 'LoaDefnListbox'
        tag = 'Loa';
        str{1} = handles.LineCandidates;
    case 'MarmDefnListbox'
        tag = 'Marm';
        str{1} = get(handles.MarmObj1Popup,'String');
        str{2} = get(handles.MarmObj2Popup,'String');
end

% Create first popup:
val1 = find(strcmpi(c{2},str{1}));
hp1 = uicontrol(hf,'Style','popupmenu',...
    'Tag','Popup1',...
    'Position',[x y w h],...
    'String',str{1},...
    'Value', val1);

% Create second popup:
x = x+w+sep/2;
if numel(str)==2  % For: {'HaxDefnListbox', 'MarmDefnListbox'}
    val2 = find(strcmpi(c{3},str{2}));
    hp2 = uicontrol(hf,'Style','popupmenu',...
        'Tag','Popup2',...
        'Position',[x y w h],...
        'String',str{2},...
        'Value', val2);
    x = x+w+sep;
    w = tw;
else    % For: 'LoaDefnListbox'
    % Don't add extra popup
    w = bw;
end

% Add OK button
hb = uicontrol(hf,'Style','pushbutton',...
    'Position',[x y w h],...
    'Tag','OkButton',...
    'String','OK');

% Resize & position figure:
fw = x+w+1.5*sep;
set(hf,'Position',[300,200,fw,55])
positionOver(hf,gcbf,'center')

% Configure behaviour
set(hf,'CloseRequestFcn',@(a,b)delete(gcbf))
set(hb,'Callback',@(varargin)set(gcbf,'Visible','off'));

% Show figure
set(hf,'Visible','on');

% Wait for finish:
waitfor(hf,'Visible','off')

% Silent exit on window X:
if ~ishandle(hf)
    return
end

% Update from return values:
new1 = str{1}{get(hp1,'Value')};
S{sel} = strrep(S{sel},c{2},new1);
if numel(str)==2
    new2 = str{2}{get(hp2,'Value')};
    S{sel} = strrep(S{sel},c{3},new2);
end
set(hobj,'String',S)

% Close figure
delete(hf)        


% ------------------------------------------------------------------------
function cstr = editList(opt,list,cstr,varargin)
% Edit one of the following lists to keep it updated:
%
%   'Hax'   Helical Axis
%   'LoA'   Line of Action 
%   'Marm'  Moment Arm
%
% Usage:
%   
%   cstr = editList('add',list,cstr,tag)
%   cstr = editList('add',list,cstr,tag1,tag2)
%       Builds a list item from TAG or TAG1 & TAG2 (depending on LIST) and
%       appends it to the cell list in CSTR which is the contents of the
%       list LIST.  LIST is a string from the abbreviations above.
%
%   cstr = editList('remove',list,cstr,j)
%       Remove the J-th item in the cell list in CSTR and re-enumerate the
%       item tags to keep them strictly monotonically increasing
%

switch opt
    case 'add'
        % Add another item to the specified list:
        obj1 = varargin{1};         % Get first item
        if numel(varargin) == 2     % If second item exists
            obj2 = varargin{2};     % Get it also
        end
        
        % Now create the list string
        k = size(cstr,1)+1; % Item number / list position
        sep = '&bull;';     % separator for two-item definitions
        switch lower(list)
            case 'hax'
                str = sprintf('<html><b>Axis_%d:</b> %s %s %s</html>',k,obj1,sep,obj2);
                
            case 'loa'
                str = sprintf('<html><b>Line_%d:</b> %s</html>',k,obj1);
                
            case 'marm'
                str = sprintf('<html><b>Arm_%d:</b> %s %s %s</html>',k,obj1,sep,obj2);
                
        end
        % Now append it to the end
        cstr{k,1} = str;
        
    case 'remove'
        % Remove the specified list item
        j = varargin{1};
        cstr(j) = [];
        cstr = re_enumerate(cstr);
end


% ------------------------------------------------------------------------
function cstr2 = re_enumerate(cstr)
%RE_ENUMERATE Re-numer a cell list of items

if isempty(cstr)
    cstr2 = cstr;
    return
end

% First rip the label:
[tok,remain] = strtok(cstr,':');
base = strtok(tok,'_');

% Re-enumerate the labels:
for j = 1:numel(base)
    base{j} = [base{j}, '_', num2str(j)];
end

% Re-form the string:
cstr2 = strcat(base',remain')';


% ------------------------------------------------------------------------
function c = html2cell(str)
%HTML2CELL Strip html from a string and truncate into its components:

c = {};  

% Strip html:
str = strrep(str,'<html>','');
str = strrep(str,'</html>','');
str = strrep(str,'<b>','');
str = strrep(str,'</b>','');

% Get label
[tag,rem1] = strtok(str,':');

% Get first item:
[item1,rem2] = strtok(rem1,'&');
item1 = strtok(item1,':');          % Trim leading ":"

% Get Second item:
[sep,item2] = strtok(rem2,';'); %#ok<ASGLU>
item2 = strtok(item2,';');          % Trim leading ";"

% Build output from components that exist.

% Tag:
c{1} = strtrim(tag);

% Item 1:
c{2} = strtrim(item1);

% Item 2 (if it exists)
if ~isempty(item2) || ( iscell(item2) && ~all(cellfun(@isempty,item2)) )
    c{3} = strtrim(item2);
end


% ------------------------------------------------------------------------
function str = current_item(hObj)
% Function to easily get the string of the current list item
cstr = get(hObj,'string');
str = cstr{get(hObj,'Value')};


% ------------------------------------------------------------------------
function html_str = ttgen(opt)
%TTGEN Generate HTML tooltip string
%
%   We do this in code instead of through GUIDE because it's easier to
%   format long multi-line strings here.  
%

switch upper(opt)
    case 'HAX'
        html_str = ['<html><b>Helical Axis</b><br>'...
            'A helical axis is defined by the motion one object<br>'...
            'relative to another.',...
            '</html>'];
        
    case 'LOA'
        html_str = ['<html><b>Line of Action</b><br>',...
            'A line of action defines the axis of a force which acts<br>',...
            'to produce rotation about the helical axis.',...
            '</html>'];
        
    case 'MARM'
        html_str = ['<html><b>Moment Arm</b><br>',...
            'The moment arm of a force is the shortest distance between<br>',...
            'the line of action of the force and the axis about which<br>',...
            'the object is rotating.  Additionally, the moment arm is<br>'...
            'also mutually perpendicular to the two vectors.',...
            '</html>'];
        
end


% ------------------------------------------------------------------------
function updateHaxPopups(schemaORhobj,eventdata)
%Keep Helical Axis popup menus in line with the components of the currently
%seleced item in the Helical Axis definition listbox
%
% Manual call, either of:
%   updateHaxPopups(handles.HaxDefnListbox)
%
% Listener call:
%   updateHaxPopups(schema,eventdata)

if nargin == 1              % Manual call
    hobj = schemaORhobj;
else                        % Listener call
    hobj = double(eventdata.AffectedObject);
end

if isempty(get(hobj,'String'))  % Nothing in the list
    return
end

handles = guidata(hobj);

% Get selection & contents:
val = get(hobj,'Value');
C = html2cell(get(hobj,'String'));

% Configure HaxObj1Popup:
c1 = get(handles.HaxObj1Popup,'String');
id1 = find(strcmpi(c1,C{2}(val)));
set(handles.HaxObj1Popup,'Value',id1);

% Configure HaxObj2Popup:
c2 = get(handles.HaxObj2Popup,'String');
id2 = find(strcmpi(c2,C{3}(val)));
set(handles.HaxObj2Popup,'Value',id2);


% ------------------------------------------------------------------------
function updateMarmPopups(schemaORhobj,eventdata)
% Manual call, either of:
%   updateMarmPopup(handles.HaxDefnListbox)
%   updateMarmPopup(handles.LoaDefnListbox)
%
% Listener call:
%   updateMarmPopup(schema,eventdata)

if nargin == 1              % Manual call
    hobj = schemaORhobj;
else                        % Listener call
    hobj = double(eventdata.AffectedObject);
end
handles = guidata(hobj);

% Get/configure MarmDefnListbox String:
defn_cell = html2cell(get(hobj,'String'));
value     = get(hobj,'Value');
switch get(hobj,'Tag')
    case 'HaxDefnListbox'
        target = handles.MarmObj1Popup;
    case 'LoaDefnListbox'
        target = handles.MarmObj2Popup;
    otherwise
        error('Unsupported')
end

% Check/set string:
contents = defn_cell{1};
if isempty(contents)
    contents = '--';   % Can't be empty
end
set(target,'String',contents,'Value',value)

% ------------------------------------------------------------------------
function MarmListUpdater(handles,oldNames,newNames)
%MARMLISTUPDATER Keep the moment arm definition list up to date
%
%   The main problem with keeping the Moment arm list up to date is that
%   when the user deletes a helical axis or a line of action, any higher
%   items in those lists get re-numbered.  So here we take care of both
%   removing any moment arms that reference deleted definitions, and also
%   replacing the higher elements with their new defintions.
%
% Note: because this function uses a string replace action, items may only
% be enumerated in the range 0-9, otherwise tags could get re-labelled
% incorrectly.  Might need to make this more robust if we want more than 9
% definitions (which is not likely).

% Get the old Moment arm defnitions:
mcstr = get(handles.MarmDefnListbox,'String');

if isempty(mcstr)
    return
end

% If NEWNAMES has a missing item, then any line in the cell array
% MCSTR that contains the corresponding old item must be removed
for j = 1:numel(oldNames)
    
    % If NEWNAMES has a missing item, then we need to remove any MCSTR
    % lines that contain the corresponding old item
    if isempty(newNames{j})
        k = regexp(mcstr,oldNames{j});      % Find the old name
        mcstr(~cellfun(@isempty,k)) = [];   % Remove any lines containing it
    else
        mcstr = strrep(mcstr,oldNames{j},newNames{j}); 
    end    
    
end

% Fix up the numbering:
mcstr = re_enumerate(mcstr);

% Now publish to the listbox:
set(handles.MarmDefnListbox,'String',mcstr)


% --- Executes on selection change in MarmObj1Popup.
function MarmObj1Popup_Callback(hObject, eventdata, handles)
% hObject    handle to MarmObj1Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MarmObj1Popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MarmObj1Popup


% --- Executes during object creation, after setting all properties.
function MarmObj1Popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MarmObj1Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MarmObj1Popup.
function MarmObj2Popup_Callback(hObject, eventdata, handles)
% hObject    handle to MarmObj1Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MarmObj1Popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MarmObj1Popup


% --- Executes during object creation, after setting all properties.
function MarmObj2Popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MarmObj1Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function varargout = LoadDicomUtil(varargin)
% LOADDICOMUTIL MATLAB code for LoadDicomUtil.fig
%
%   USAGE - blocking behaviour:
%       LoadDicomUtil
%       files = LoadDicomUtil
%       files = LoadDicomUtil({seed_path})
%
%   USAGE - non-blocking behaviour;
%       hfig = LoadDicomUtil({},'-detach')
%       hfig = LoadDicomUtil({seed_path},'-detach')
%       

% Last Modified by GUIDE v2.5 15-Nov-2012 13:18:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LoadDicomUtil_OpeningFcn, ...
                   'gui_OutputFcn',  @LoadDicomUtil_OutputFcn, ...
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
function LoadDicomUtil_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before LoadDicomUtil is made visible.
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LoadDicomUtil (see VARARGIN)

% Ensure dependendencies are included:
addpath(genpath('lib/2d'))
addpath(genpath('lib/shared'))

% Get input arguments
%   Note that the first argument in varargin{} can't be a string, because
%   the main gui function (above) will try to inrpret it as an attempt to
%   run a subfunction.  So the first input (default path) should be wrapped
%   as a cell.  Perhaps there is a better way to deal with this...?
pth = pwd;
if numel(varargin) > 0  && ~isempty(varargin{1}) && ~isempty(varargin{1}{1})
    pth = varargin{1}{1};   % Path should be wrapped as a cell.
end

% Position gui over the top of calling gui
positionOver(handles.figure1,gcbf);

% Configure axes:
set(handles.axes1,'Visible','off')          % Don't display axes
imshow([],'Parent',handles.axes1)           % Show nothing
set(handles.Slider_Frame, 'Visible','off')  % Don't display slider


addlistener(handles.Listbox_Files,'String','PostSet',@list_listener);
addlistener(handles.Listbox_Files,'Value','PostSet',@list_listener);

set_directory(pth)

%set([handles.Pushbutton_OpenFolder, handles.Pushbutton_UpOneDir],...
%    'BackgroundColor',[0.95 0.95 0.95])
% Set some icons:
load icons;
set(handles.Pushbutton_OpenFolder,'CData',icons.explore_dir)
set(handles.Pushbutton_UpOneDir,'CData',icons.up_dir)
set(handles.Pushbutton_RefreshDir,'CData',icons.refresh)
%flatten_buttons([...
%   handles.Pushbutton_OpenFolder,...
%   handles.Pushbutton_UpOneDir,...
%   handles.Pushbutton_RefreshDir]);

% Default command line output
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Now manage the calling method & command line outputs:
%   1) if called with flag '-detach', return handle to the command line, and
%       outputs should be obtained using a listner (see
%       closeRequestFunction)
%   2) if flag not present, wait for closure and return result to
%       command line
if any( strcmpi(varargin,'-detach') ) % Check for -detach flag
	% Default command line output already set (above)
else
    % See OutputFcn & CloseReqestFcn before modifying these lines
    set(handles.figure1,'WindowStyle','modal')  % Blocking; also referenced in OutputFcn
    uiwait(handles.figure1);                    % Wait for uiresume before returning
end

end


% ------------------------------------------------------------------------
function varargout = LoadDicomUtil_OutputFcn(hObject, eventdata, handles) 
% --- Outputs from this function are returned to the command line.
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% If it was a modal & blocked by uiwait, then we have arrived here by
% uiresume and we can delete the figure and exit:
if strcmpi( get(handles.figure1,'WindowStyle'), 'modal' ) 
    delete(handles.figure1)                 % Delete figure, return to command line
end

end


% ------------------------------------------------------------------------
function Listbox_Files_Callback(hObject, eventdata, handles)
% --- Executes on selection change in Listbox_Files.
% hObject    handle to Listbox_Files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Listbox_Files contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Listbox_Files

contents = cellstr(get(hObject,'String'));
value = contents{get(hObject,'Value')};

% Drop pre-pended space or '+'
if ~any( strcmpi(value,{'.','..'}) )
    value = value(2:end);   
end

cwd = get_directory();

OPEN  = isequal('open', get(handles.figure1,'SelectionType'));
ISDIR = exist([cwd value],'dir')==7;

switch value
    case '.'    % Current directory
        % do nothing
    case '..'   % Parent directory
        % move up a directory
        if OPEN
            set_directory( parent_dir_of(cwd) )
        end 
        
    otherwise   % File or directory:
        
        if ISDIR &&  OPEN % Open folder    
            set_directory([cwd value filesep])
            
        end
        % Listeners handle everything else
        
        %{
        elseif ISDIR     % Select folder for preview
            %disp('Select folder')
            %dcm_preview([cwd value filesep],handles.axes1)
            
        elseif ~ISDIR    % Select file for preview
            %disp('Select file')
            %dcm_preview([cwd value],handles.axes1)
            
        end
        %}
        
        
end

guidata(hObject,handles)

end


% ------------------------------------------------------------------------
function Listbox_Files_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to Listbox_Files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% ------------------------------------------------------------------------
function Checkbox_CacheDICOMDIR_Callback(hObject, eventdata, handles)
% hObject    handle to Checkbox_CacheDICOMDIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end


% ------------------------------------------------------------------------
function Checkbox_FilterDICOMDIR_Callback(hObject, eventdata, handles)
% --- Executes on button press in Checkbox_FilterDICOMDIR.
% hObject    handle to Checkbox_FilterDICOMDIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filter_checkbox_handler(hObject)

end


% ------------------------------------------------------------------------
function Checkbox_FilterExtension_Callback(hObject, eventdata, handles)
% --- Executes on button press in Checkbox_FilterExtension.
% hObject    handle to Checkbox_FilterExtension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filter_checkbox_handler(hObject)

end


% ------------------------------------------------------------------------
function Edit_Extensions_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_Extensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edit_Extensions as text
%        str2double(get(hObject,'String')) returns contents of Edit_Extensions as a double

% On return, refresh file list:
file_list_refresh(get_directory(handles))

end


% ------------------------------------------------------------------------
function Edit_Extensions_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to Edit_Extensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% ------------------------------------------------------------------------
function Listbox_Info_Callback(hObject, eventdata, handles)
% --- Executes on selection change in Listbox_Info.
% hObject    handle to Listbox_Info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Listbox_Info contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Listbox_Info

end


% ------------------------------------------------------------------------
function Listbox_Info_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to Listbox_Info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'))

end


% ------------------------------------------------------------------------
function Popup_Path_Callback(hObject, eventdata, handles)
% --- Executes on selection change in Popup_Path.
% hObject    handle to Popup_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Popup_Path contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Popup_Path

% Get selected item
contents = cellstr(get(hObject,'String'));
sel = get(hObject,'Value');

% Build path string:
c = flipud(contents(sel:end))';
c(2,:) = {filesep};
s = [c{:}];
if ~ispc        % on unix systems
    s(1) = [];  % drop the additional preceding '/'
end

% Refresh display:
set_directory(s)


end


% ------------------------------------------------------------------------
function Popup_Path_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to Popup_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% ------------------------------------------------------------------------
function Pushbutton_Browse_Callback(hObject, eventdata, handles)
% --- Executes on button press in Pushbutton_Browse.
% hObject    handle to Pushbutton_Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

p = get_directory();
dname = uigetdir(p,'Select directory');
if isequal(dname,0)
    return
end
set_directory(dname)

end


% ------------------------------------------------------------------------
function Pushbutton_OpenFolder_Callback(hObject, eventdata, handles)
% --- Executes on button press in Pushbutton_OpenFolder.
% hObject    handle to Pushbutton_OpenFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ismac
    system(['open "' get_directory() '"'] );
elseif ispc
    system(['explorer "' get_directory() '"']);
else
    try
        system(['nautilus "' get_directory() '"']);
    catch ME
        disp('Not available on linux.  Error reported:')
        disp(ME)
    end        
end

end


% ------------------------------------------------------------------------
function Pushbutton_RefreshDir_Callback(hObject, eventdata, handles)
% --- Executes on button press in Pushbutton_RefreshDir.
% hObject    handle to Pushbutton_RefreshDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_directory( get_directory() )
        
end


% ------------------------------------------------------------------------
function Pushbutton_UpOneDir_Callback(hObject, eventdata, handles)
% --- Executes on button press in Pushbutton_UpOneDir.
% hObject    handle to Pushbutton_UpOneDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_directory( parent_dir_of( get_directory() ) )

end


% ------------------------------------------------------------------------
function Pushbutton_LoadFiles_Callback(hObject, eventdata, handles)
% --- Executes on button press in Pushbutton_LoadFiles.
% hObject    handle to Pushbutton_LoadFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Everything happens in here:
figure1_CloseRequestFcn(hObject, [], guidata(hObject) )

end


% ------------------------------------------------------------------------
function Pushbutton_Cancel_Callback(hObject, eventdata, handles)
% --- Executes on button press in Pushbutton_Cancel.
% hObject    handle to Pushbutton_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Everything happens in here:
figure1_CloseRequestFcn(hObject, [], guidata(hObject) )

end 


% ------------------------------------------------------------------------
function Slider_Frame_Callback(hObject, eventdata, handles)
% hObject    handle to Slider_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


imglist = study_images_from_list(handles.Listbox_Info);

val = get(hObject,'Value');
k = round(val);

update_preview(imglist{k},handles.axes1)

end

% ------------------------------------------------------------------------
function Slider_Frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slider_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end

% ------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% --- Executes when user attempts to close figure1.
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Different closing methods / outputs are determined by how the figure is
% closed:

switch hObject
    
    case handles.figure1
        % In this case, user closed using the window X button
        % OUTPUT: []
        handles.output = [];   
        
    case handles.Pushbutton_Cancel
        % User pressed cancel.
        % OUTPUT: []
        handles.output = [];    
        
    case handles.Pushbutton_LoadFiles
        
        % Get list of files to load - abort when no valid file selected
        f = get_current_dicom_files(handles);
        if isempty(f)
            warndlg('No valid files selected!','No files','modal')
            return
        end
        handles.output = f;
        
    otherwise
        return
end


% Update handles so listener gets up-to-date info:
guidata(hObject,handles)

% Now switch the two methods - using uiwait, or not:
dstack = dbstack;
if any(strcmpi('uiwait',{dstack.name}))
    % closing with the output sent to the command line
    % OutputFcn will delete the figure for us
    uiresume(handles.figure1)

else
    % Closing with the output being grabbed from handles.output by a
    % listener
    
    % First, hide the figure so the user can't mess anything up:
    set(handles.figure1,'Visible','off')
    
    % Next, trigger the listener from invoking program, MRIMagic:
    set(handles.figure1,'HitTest','off')
    
    % Then when the listener completes, it will return here
    delete(handles.figure1)
end

end %figure1_CloseRequestFcn()



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%      GUI REFRESH / UPDATE FUNCTIONS
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------------------
function file_list_refresh(pathstring,handles)
% Update Listbox_Files from specified directory
% Usage:
%   file_list_refresh(pathstring)
%   file_list_refresh(pathstring,handles)

if exist('handles','var') == 0
    handles = gethandles();
end

hObject = handles.Listbox_Files;    % This is the object to operate on
D = dir(pathstring); 
if isempty(D)
    errordlg('This folder is inaccessible','Directory error')
    return
end

% Sort directories first:
d = D([D.isdir]);
dlist = {d.name}';

% Apply filters to file list:
if get(handles.Checkbox_FilterDICOMDIR,'Value')
    f = dir([pathstring 'DICOMDIR']);
elseif get(handles.Checkbox_FilterExtension,'Value')
    ex = get(handles.Edit_Extensions,'String');
    ex = extStr2extCell(ex);
    f = cellfun(@(e)dir([pathstring e]),ex,'UniformOutput',false);
    f = cat(1,f{:});
else 
    % Otherwise include all files:
    f = D(~[D.isdir]);
end
flist = {f.name}';

% Don't show hidden files
flist = flist(~cellfun(@(p)any(p==1),strfind(flist,'.')));

% Don't show hidden directories, but do show . and ..
dhidden = cellfun(@(p)any(p==1),strfind(dlist,'.'));
dhidden(1:2) = false;   % These are . and ..
dlist = dlist(~dhidden);

% Prepend space or plus:
if numel(dlist) > 2
    dlist(3:end) = cellsmash('+',dlist(3:end));
end
flist = cellsmash(' ',flist);

% Concatenate directories & files:
list = [dlist;flist];

% Now update - in the correct order to prevent warnings
set(hObject,'ListboxTop', min( [get(hObject,'ListboxTop'), numel(list)] ) )
drawnow update
set(hObject,...
    'String', list,...
    'Value',      min( [get(hObject,'Value'),numel(list)] ) )
end %file_list_refresh()


% ------------------------------------------------------------------------
function filter_checkbox_handler(hObject)
% Callback for managing filter checkboxes
% See also FILTER_FILE_LIST()
handles = guidata(hObject);
hfext = handles.Checkbox_FilterExtension;
hddir = handles.Checkbox_FilterDICOMDIR;

% Manage mutual exclusivity
if (hObject == hfext) && (get(hObject,'Value')==1)
    set(hddir,'Value',0)
elseif (hObject == hddir) && (get(hObject,'Value')==1)
    set(hfext,'Value',0)
end

% Refresh file list with new settings:
file_list_refresh(get_directory(handles))

end %filter_checkbox_handler()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%      LISTENER FUNCTIONS
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------------------
function list_listener(~,ev)
% This function needs to be called anytime the currently selected file
% changes.  This happens under the two events:
%   - Selection change
%   - String (list contents) update
%
% But selection is linked to the size of the list, so we get duplicate
% callbacks when contents are changed and the selection is then moved to be
% within the length of the list.  We catch this by testing to see if
% selection is outside the list length, if it is, there will be an ensuing
% call to this listner, so we can skip the present one.

%prop = get(sch,'Name') % Property that triggered the listener { 'String' | 'Value' }
obj = ev.AffectedObject;
contents = cellstr(get(obj,'String'));

% Test to see if selection is out-of-bounds.  If so, skip this call and act
% on the ensuing call
if get(obj,'Value') > numel(contents)
    % We can expect a value change, so skip this call
    return
end

% Ok, now we shouldn't get double-calls here.

handles = guidata(obj);

% Now farm out the operation:
obj = double(obj);  % Convert to its double counterpart, which is the norm
switch obj
    
    case handles.Listbox_Files
        % Get current path/file:
        fname = fullpath_of_selection(handles);
        
        % If it's a DICOMDIR, do this:
        %   - make study list
        %   - update preview
        
        % If it's a DICOM image file, do this
        %   - make header list
        %   - update preview
        
        update_infolist(fname, handles.Listbox_Info, handles)
        
        % Try to display it:
        update_preview(fname, handles.axes1);
        
        
    otherwise
        error('unhandled')
        
end

end %list_listener()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%      SUPPORT FUNCTIONS
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------------------
function C = cellsmash(A,B)
% Smash cell arrays A & B into the same cell array
% Both A & B can be cell arrays of the same size, or either A or B can be a
% cell array with the other being a string for insertion into every element
% Usage:
%   C = cellsmash({'hi';'bye';'fly'}, ' Barry')
%   C = cellsmash('Barry ',{'thought','said','did'})
%   C = cellsmash({'Barry';'Jo';'Bill'},{' can'; ' did'; ' will'})
%   C = cellsmash({'hi ','Barry';'bye ','Jo'})
if nargin == 1  % cellsmash(A)
    for j = 1:size(A,1)
        C{j,1} = [A{j,:}];
    end
elseif iscell(A) && iscell(B)
    assert(isequal(numel(A),numel(B)),'If two cell arrays are provided, they must be the same size')
    C = cellfun(@(a,b)[a,b],A,B,'UniformOutput',false);
elseif iscell(A) && ~iscell(B)
    C = cellfun(@(a) [a,B],A,'UniformOutput',false);
elseif ~iscell(A) && iscell(B)
    C = cellfun(@(b) [A,b],B,'UniformOutput',false);
else
    error('Bad inputs to cellsmash()')
end


end %cellsmash()


% ------------------------------------------------------------------------
function exlist = extStr2extCell(exstr)
% Convert user-defined extension string to cell list
d = ';';
ex = lower(exstr);
ex = regexprep(ex,',',d); % Convert commas to delimiters
ex = regexprep(ex,' ',d); % Convert spaces to delimiters
excell = textscan(ex,'%s','Delimiter',d,'MultipleDelimsAsOne',true);
exlist = excell{1};

end %extStr2extCell()


% ------------------------------------------------------------------------
function flatten_buttons(hobjs)
% This makes prettier buttons by removing their border.  But there are
% drawbacks:
%   - The figure must first be visible (which is forced by this function)
%   - Setting the border with java seems to be painfully slow
if exist('findjobj','file')~=2
    return
end
% Buttons need to be drawn before they can be modified with java...
set(ancestor(hobjs(1),'figure'),'Visible','on')
drawnow
for j = 1:numel(hobjs)
    jh = findjobj( hobjs(j) );
    jh.setBorder([])
end

end %flatten_buttons()


% ------------------------------------------------------------------------
function handles = gethandles()
% Get the guidata without having to provide the figure handle.
% A convenience thing so we don't have to pass handles around everywhere,
% and so we always get the current handles data.
persistent hf
if isempty(hf) || ~ishandle(hf)
    hfigs = findall(0,'Type','Figure');
    fnames = get(hfigs,'FileName');
    if iscell(fnames)
        hf = hfigs( ~cellfun(@isempty,regexp(fnames,mfilename)) );
    else
        hf = hfigs; % Only one figure; this must be it
    end
end
handles = guidata(hf);

end %gethandles()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%      FILE / PATH HELPER FUNCTIONS
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------------------
function fname = filename_of_selection(handles)
if exist('handles','var') == 0
    handles = gethandles();
end
obj = handles.Listbox_Files;
contents = cellstr(get(obj,'String'));
value = contents{get(obj,'Value')};
fname = value;
% Now chop the first character off, because it is either a space or a '+':
if any(strcmpi(value(1),{'+',' ',filesep}))
    fname = value(2:end);
end

end %filename_of_selection()


% ------------------------------------------------------------------------
function fullname = fullpath_of_selection(handles)
% Get full path/filename of currently selected item in Listbox_Files
if exist('handles','var') == 0
    handles = gethandles();
end
fullname = [get_directory(handles), filename_of_selection(handles)];

end %fullpath_of_selection()


% ------------------------------------------------------------------------
function f = get_current_dicom_files(handles)
% Get curent files & return in cell array:
obj = handles.Listbox_Info;
ud = get(obj,'UserData');

% F will only be populated by DICOMDir objects, or DICOM image file paths
if isempty(ud)
    f = {};
    
elseif isa(ud,'DICOMDir')
    f = study_images_from_list(obj);    % Get list of files in study
    
elseif isa(ud,'char') && isdicom(ud)    
    f = {ud};                           % Get single dicom filepath
    
end

end %get_current_dicom_files()


% ------------------------------------------------------------------------
function p = get_directory(handles)
% Get path of current directory
if exist('handles','var') == 0
    handles = gethandles();
end
% It's stored as the string of this object:
p = get(handles.Text_Pathstring,'String');

end %get_directory()


% ------------------------------------------------------------------------
function parent = parent_dir_of(pathstring)
% Get the parent directory of specified directory
if isequal(pathstring(end),filesep)
    n = 2;
else
    n = 1;
end
ids = find(pathstring == filesep,n,'last');
parent = pathstring(1:ids(1));

end %parent_dir_of()


% ------------------------------------------------------------------------
function path_list_refresh(pathstring,handles)
% Usage:
%   path_list_refresh(pathstring)
%   path_list_refresh(pathstring,handles)

if exist('handles','var') == 0
    handles = gethandles();
end

% Chop into segments & reverse order:
c = [filesep; path2cell(pathstring)'];
c = flipud(c);

% Set the contents & selection value:
set(handles.Popup_Path,'String',c)
set(handles.Popup_Path,'Value',1)


end %path_list_refresh()


% ------------------------------------------------------------------------
function set_directory(pathstring)
% Central function for setting the currently displayed directory

% First check if the directory is accessible:
if isempty( dir(pathstring) ) % dir() should always return at least two entries: . and ..
    errordlg('This folder is inaccessible','Directory error')
    return
end

% Ensure it contains the trailing filesep
if ~isequal(pathstring(end),filesep)
    pathstring(end+1) = filesep;
end

handles = gethandles();

% Set the path text string
set(handles.Text_Pathstring,'String',pathstring)

% Refresh the path dropdown list:
path_list_refresh(pathstring,handles)

% Refresh the file list:
file_list_refresh(pathstring,handles)

end %set_directory()



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%      PREVIEWING FUNCTIONS & THEIR HELPERS
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
function update_preview(pth,axs)
% Preview the specified image on axis axs
I = [];
map = [];
try
    if isdir(pth)
        % do nothing for paths - display as empty
        
    elseif isdicom(pth) && ~isdicomdir(pth)
        I = dicomread(pth);
        
    elseif isimage(pth)
        [I,map] = imread(pth);
        
    end
catch ME
    fprintf('Image display failed...\n');
    fprintf('      error: %s\n',ME.identifier);
    fprintf('    message: %s\n',ME.message);
end
imshow(I,map,'Parent',axs)


end %update_preview()


% ------------------------------------------------------------------------
function update_infolist(pth,hObject,handles)
% This function updates the info list, handles.Listbox_Info, which is
% hObject
MONOSPACED = 'Monospaced';
NORMAL = 'Helvetica';
font = MONOSPACED;      % Set default
cbk = [];
ud  = [];
S = [];
v = get(hObject,'Value');
lbt = get(hObject,'ListBoxTop');
slideVis = 'off';


if isdir(pth)
	% Normal directories, do nothing; but DICOM directories with a DICOMDIR
	% file present, handle them the same as 
    
elseif isdicomdir(pth)
    usecache = get(handles.Checkbox_CacheDICOMDIR,'Value');
    [S,ud] = list_from_dicomdir(pth,usecache);
    cbk = @study_preview;
    uicontrol(hObject)      % switch focus to the study list
    lbt = 1;
    slideVis = 'on';
    
elseif isdicom(pth)
    font = NORMAL;
    S = header_from_dicom(pth,font);
    ud = pth;
    if lbt > size(S,1)
        lbt = size(S,1);
    end
elseif isimage(pth)
    S = header_from_image(pth,font);
    
    % do nothing - show blank
end
if isempty(v) || v == 0
    v = 1;
elseif v > size(S,1)
    v = size(S,1);
end
set(hObject,...
    'Callback', cbk,...
    'String',   S,...
    'FontName', font,...
    'ListboxTop', lbt,... 
    'Value',    v,...
    'UserData', ud)

% Hide/show slider:
set(handles.Slider_Frame,'Visible',slideVis)


end %update_infolist()


% ------------------------------------------------------------------------
function study_preview(hObject,~)

% hObject is now the handle to the listbox
imglist = study_images_from_list(hObject);

% Start by selecting the middle image:
n = numel(imglist);
k = (n/2);
s = [1 10]/n; 

handles = guidata(hObject);

% Set the slider:
set(handles.Slider_Frame,'Min',1,'Max',n,'Value',k,'SliderStep',s)

% Display the image:
k = round(k);
update_preview(imglist{k}, handles.axes1);

end %study_preview()


% ------------------------------------------------------------------------
function tf = isdicomdir(pname)
% Test if file is a dicomdir file
[~,f,~] = fileparts(pname);
% Use IF rather than direct logic to make short-circuiting clearer:
if exist(pname,'file')==2 && ...   % Check file exists
        isdicom(pname) && ...            % Should pass this too
        isequal(f,'DICOMDIR')          % Name should be DICOMDIR
    tf = true;
else
    tf = false;
end
end %isdicomdir()


% ------------------------------------------------------------------------
function tf = isdir(pname)
tf = exist(pname,'dir')==7;
end %isdir()


% ------------------------------------------------------------------------
function tf = isimage(fname)
% Test to see if file is probably an image (by examining file extension)
[~,~,ex] = fileparts(fname);
formats = imformats;
imexts = {formats.ext};
imexts = [imexts{:}];
imexts = cellsmash('.',imexts);
tf = any( strcmpi(imexts, ex) );
end %isimage()


% ------------------------------------------------------------------------
function H = header_from_dicom(dcmimgfile,font)
props = {...
    'RequestedProcedureDescription',... % OR 'PerformedProcedureStepDescription'
    'RequestedProcedureComments',...
    'AcquisitionDate',...
    'PatientID','PatientName',...
    'StudyDescription','SeriesDescription','SeriesNumber',...
    'BodyPartExamined','Modality',...
    'AcquisitionNumber','InstanceNumber',...
    'RepitionTime','EchoTime','NumberOfAverages','ImagingFrequency',...
    'NumberOfSlices','SliceThickness','SliceLocation','PixelSpacing',...
    'ImageType','Filename','Rows','Columns'...
    };
props = sort(props);
i = dicominfo(dcmimgfile);
% Copy all requested fields that are present into cell array of chars:
C = {};
sep = ':  ';
for j = 1:numel(props)
    if isfield(i,props{j})
        p = props{j};
        v = header_data_sanitize(i.(props{j}));
        C(end+1,:) = { p , sep, v } ; 
        %h.(props{j}) = i.(props{j});
    end
end

H = cell_html_formatted_list(C,font);

end %header_from_dicom()

% ------------------------------------------------------------------------
function L = cell_html_formatted_list(C,font)
% Find the longest item in each column, & calculate the width (in pixels) of each
% column:
[nr,nc] = size(C);
cwidths = NaN(1,nc);
hf = figure('Visible','off');
for k = 1:nc
    hu = uicontrol(hf,...
        'Style','text',...
        'HorizontalAlignment','Left',...
        'FontName',font,...
        'String',C(:,k));
    extent = get(hu,'Extent');
    cwidths(k) = extent(3);
end
close(hf);
for k = 1:nr
    L{k,1} = html_columnised_row(cwidths,C(k,:));
end

    % -----------------------------------------
    function htmlrow = html_columnised_row(widths,entries)
        nc = numel(widths);
        pad = 2; % [pixels]
        widths = widths+pad;
        table_width =  sum(widths);
        pre1 = sprintf('<html><table width=%s style=''table-layout:fixed''>', num2str(table_width));
        pre2 = [];
        content = '<tr>';
        for j = 1:nc
            pre2 = [pre2 sprintf('<col width=%s>', num2str(widths(j)))];  %#ok<*AGROW>
            content = [content sprintf('<td width="%s">%s</td>',num2str(widths(j)),entries{j})]; 
        end
        content = [content '</tr>'];
        post = '</table></html>';
        htmlrow = [pre1 pre2 content post];
    end %html_columnised_row()

end %cell_html_formatted_list()


% ------------------------------------------------------------------------
function s = header_struct2str(S)
fields = fieldnames(S);
s = '';
for k = 1:numel(fields)
    fkname = fields{k};
    fkdata = S.(fkname);
    if isa(fkdata,'struct')
        fkdata = header_struct2str(fkdata);
    elseif isnumeric(fkdata)
        fkdata = num2str(fkdata);
    end
    if k == 1
        fsep = '';
    else
        fsep = ',  ';  % field separator
    end
    s = [s fsep fkname ': ' fkdata];
end
end %header_struct2str()

% ------------------------------------------------------------------------
function H = header_from_image(impath,font)
% Create header using image information gathered using IMINFO
% See also header_from_dicom()

[pth,name,ext] = fileparts(impath);

% Get info & structure header text:
info = imfinfo(impath);

% Make some adjustments:
info.Filename = [name ext];
info.Path = pth;
info.Format = upper(info.Format);

% Create list of properties we want to display:
props = {'Filename','Path','Format','Width','Height',...
    'BitDepth','ColorType','Compression'};

% Create cell array of relevant data:
C = {};
sep = ':  ';
for p = props
    if isfield(info,p)
        % See also header_from_dicom
        v = header_data_sanitize(info.(p{1}));        
        C(end+1,:) = { p{1}, sep, v };
    end
end

H = cell_html_formatted_list(C,font);

end %header_from_image()

% ------------------------------------------------------------------------
function data = header_data_sanitize(data)
% Convert numeric or structure into string.  Strings pass through
% un-touched
if isnumeric(data)
    data = num2str(data(:)');     % force 1-by-n format, then convert
elseif isa(data,'struct')
    data = header_struct2str(data);
end
end %header_data_sanitize()

% ------------------------------------------------------------------------
function [L,D] = list_from_dicomdir(dcmdirfile,usecache)

handles = gethandles();

% Handle caching:
pth = fileparts(dcmdirfile);
cachefile = [pth filesep 'dicomdir.mat'];

if usecache && exist(cachefile,'file')==2
    % Load from cache - provides D
    load(cachefile)
else
    % Lock figure
    hlock = FigLocker.Lock(handles.figure1);
    hlock.addprogressbar;
    hlock.settext('Parsing DICOMDIR file.  This may take a moment...');
    % Read from DICOMDIR file
    D = DICOMDir(dcmdirfile);
    % Then save if required:
    if usecache
        save(cachefile,'D')
    end
    % Unlock
    hlock.unlock;
end


% Drop empty series: (where SeriesNumber==0)
valid_series = ~cellfun(@(n)isempty(n)|isequal(n,0),{D.series.SeriesNumber})';
D.series = D.series(valid_series);
D.seriesMap = D.seriesMap(valid_series,:);

% Find number of images in each study:
for j = size(D.seriesMap,1) : -1 : 1
    nimages(j) = numel(D.imagesInSeries(D.seriesMap(j,:)));
end

% Build cell containing all the metadata to display in list:
% (Robust with if/else in case fields don't exist)
C = cell(numel(D.series),1);
C(:,1) = {D.patients.PatientName.FamilyName};
if isfield(D.series,'Modality')
    C(:,end+1) = {D.series.Modality}';
end
if isfield(D.studies,'StudyDescription')
    C(:,end+1) = {D.studies.StudyDescription};
end
if isfield(D.series,'ProtocolName')             % PHILIPS
    C(:,end+1) = {D.series.ProtocolName};
elseif isfield(D.series,'SeriesDescription')    % SIEMENS (& others?)
    C(:,end+1) = {D.series.SeriesDescription};
end
C(:,end+1) = cellsmash('Images: ', cellstr(num2str(nimages(:))) );

% Convert to a string suitable for the listbox:
L = cell_to_listbox_string(C);


    % ------------------------------------
    function L = cell_to_listbox_string(C)
        % Re-format a N-by-M cell to a N-by-1 cell.  
        % C must contain only char elements
        [nr,nc] = size(C);
        sep = repmat(sprintf('  '),nr,1);
        for k = 1:nc
            c = char(C(:,k));
            if k == 1
                S = c;
            else
                S = [S sep c];
            end
        end
        [nr,nc] = size(S);
        L = mat2cell(S,ones(1,nr),nc);
    end %cell_to_listbox_string()

end %list_from_dicomdir() 


% ------------------------------------------------------------------------
function imgcell = study_images_from_list(hObject)
% Get the cell array of images from the currently selected study
assert(isequal(get(hObject,'Tag'),'Listbox_Info'),'This function acts on Listbox_Info')
v = get(hObject,'Value');       % Get id of selected row
D = get(hObject,'UserData');    % get DICOMDir object
assert(isa(D,'DICOMDir'))
% Get the image list using the series map:
pss = D.seriesMap(v, :);
imgcell = D.imagesInSeries( pss );

end %study_images_from_list() 



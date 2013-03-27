function save_session(hObject,handles)
% SAVE_SESSION Helper method called by both programs which handles session
% saving, including "save" and "save as"


% Lock figure
locker = FigLocker.Lock(handles.figure1);

% Get action states based on hObject:
SAVE   = isempty(regexpi(get(hObject,'Tag'),'As','ONCE'));
%SAVEAS = ~SAVE;
PATHEXISTS = ~isempty(handles.sessionPath);

% Set default path:
defaultPath = handles.sessionPath;


% Switch based on save method
if SAVE && PATHEXISTS
    % Save over the top of the specified file without prompting:
    PROMPT = false;    
else % "Save As"
    PROMPT = true;
    % Check the default path
    if isempty(defaultPath) && ~isempty(handles.userPath)
        defaultPath = handles.userPath;
    end
end

% Perform save:
[pathname,ok] = Session.save(handles,defaultPath,PROMPT);

% Update handles if save was successful:
if ok
    handles.sessionPath = pathname;
    guidata(hObject,handles)
end

% Unlock figure:
locker.unlock; 



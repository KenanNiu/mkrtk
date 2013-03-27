function [handles,ok] = load_session(handles)

ok = false;

% Lock figure:
locker = FigLocker.Lock(handles.figure1);

% Get file from usesr
[filename, pathname] = uigetfile(...
    {'*.mat','MAT-files (*.mat)'},'Load session',handles.sessionPath);

if isequal(filename,0)
    locker.unlock
    return
end

% Determine the invoking file:
ds = dbstack(1);
[~,caller,~] = fileparts(ds(1).file);

% Load the session, specifying the caller (and hence fields that should be
% loaded):
[handles,ok] = Session.load(handles,[pathname filename],caller);

% Check for success:
if ok
    % Update guidata:
    guidata(handles.figure1,handles)

end

locker.unlock;


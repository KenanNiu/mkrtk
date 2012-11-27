function runCallback(hObject,varargin)
%RUNCALLBACK Run the callback of the specified uicontrol object
%
% Usage:
%   runCallback(handles.PushButton)         eventData = [] implicit
%   runCallback(handles.PushButton, [])     eventData = [] explicit
%   runCallback(handles.PushButton, evd)    eventData = evd

cbk = get(hObject,'Callback');
switch numel(varargin)
    case 0
        cbk(hObject,[])
    case 1
        cbk(hObject,varargin{1})
    otherwise
        error('incorrect number of input arguments')
end
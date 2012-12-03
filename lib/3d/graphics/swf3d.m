function swf3d(hf,event)
%SWF3D Scroll Wheel Function for figures showing the phase-dependent 3d axes
%
% SWF3D is hung on either of the two figures, with the target uicontrol as
% shown:
%   1) Registration (main 3D program)
%       - "PhaseSlider" is a slider uicontrol which is the master control for
%         changing the displayed phase in this interface
%
%   2) Registration Solver (interface for performing the registration)
%       - "PhasePopup" is a popup uicontrol which is the master control for
%         changing the displayed phase in this interface
%
% SWF3D captures the scroll action and converts in to an incrment/decrement
% on the appropriate uicontrol.  The uicontrol's callback is then run to
% action the display refresh.


% Get scroll amount:
scroll = parseInputs(event);

% Now find the target uicontrol:
handles = guidata(hf);
if isfield(handles,'PhaseSlider')
    hObject = handles.PhaseSlider;
    valMin = get(hObject,'Min');
    valMax = get(hObject,'Max');
    
elseif isfield(handles,'PhasePopup')
    hObject = handles.PhasePopup;
    pstr = cellstr(get(hObject,'String'));
    valMin = 1;
    valMax = numel(pstr);
    
else
    return
end

%[~,fn,~] = fileparts(get(hf,'Filename'));
    
% Create function handle for main gui m-file:
%caller = str2func(fn);  

% Scroll amount & new position:
val = get(hObject,'Value');
newval = val + scroll;


if newval <= valMax && ...
        newval >= valMin

    % Set the new value
    set(hObject,'Value',newval)
    
    % Call the uicontrol's callback to action the update:
    cbk = get(hObject,'Callback');
    cbk(hObject,[])
    
else
    % Tried to scroll past the end
    return
end


% ------------------------------------------------------------------------
function scroll = parseInputs(event)
% Parse scroll amount:
%   event.VerticalScrollCount
%   event==[]
if isempty(event)
    scroll = 0;
else
    scroll = event.VerticalScrollCount;
end

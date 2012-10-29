function swf3d(hf,event)
%SWF3D Scroll Wheel Function for 3D cloud viewing
%
% SWF3D calls the slider callback function if the slider is visible

handles = guidata(hf);

%--------- Check Slider config:
% nphases = max( cellfun(@numel,handles.LoResClouds) );
% if isempty(nphases) || nphases == 0
%     nmin = 1;
%     nmax = 2;
% else
%     nmin = 1;
%     nmax = nphases;
% end
% set(handles.Slider3D,'Min',nmin);
% set(handles.Slider3D,'Max',nmax);
% 
% if isequal('off',get(handles.Slider3D,'Visible'))
%     return
% end


%--------- Adjust slider:

% Get scroll amount:
scrollAmount = parseInputs(event);

% Get the file name of the GUI (automated, in case it changes):
[~,fn,~] = fileparts(get(hf,'Filename'));
    
% Create function handle for main gui m-file:
caller = str2func(fn);  

% Scroll amount & new position:
val = get(handles.Slider3D,'Value');
newval = val + scrollAmount;


if newval <= get(handles.Slider3D,'Max') && ...
        newval >= get(handles.Slider3D,'Min')

    % Set the new value
    set(handles.Slider3D,'Value',newval)
    
    % Call the slider function in the main program:
    caller('Slider3D_Callback', handles.Slider3D, [], handles)
else
    % Tried to scroll past the end
    return
end


% ------------------------------------------------------------------------
function scrollAmount = parseInputs(event)
% Parse scroll amount:
%   event.VerticalScrollCount
%   event==[]
if isempty(event)
    scrollAmount = 0;
else
    scrollAmount = event.VerticalScrollCount;
end

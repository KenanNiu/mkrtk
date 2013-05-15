function swf2d(hf,event)
%SWF2D Scroll Wheel Function for 2D Dicom viewing
%
% SWF2D does the following:
%   - Change the currently displayed image by:
%       > Scrolling up/down to change the slice
%       > Shift + scrolling to change the phase
%   - Update necessary annotations
%
% See also TRACEWBDF, TRACEWBMF, TRACEWBUF


[hf,scrollAmount] = parseInputs(hf,event);

% Get up-to-date info:
handles = guidata(hf);
sj = current('slice',handles.axes1);
pj = current('phase',handles.axes1);

% Return if not yet configured:
if ( isempty(sj) && isempty(pj) ) || ( isempty(handles.Images) || isempty(handles.Images.X) )
    return
end

[ns,np] = handles.Images.stacksize();

% Now we check for the SHIFT modifier 
%   If shift is being held while scrolling, we are scrolling phases, if
%   not, we're scrolling slices: 
if strcmpi( 'shift', getappdata(hf,'CurrentKey') )
    pj = pj - scrollAmount;     % Increment phase
    pj = max( 1, min(np,pj) );  % Limit scrolling to indices of image stack
else
    sj = sj - scrollAmount;       % Increment scroll
    sj = max( 1, min(ns,sj) );  % Limit scrolling to indices of image stack
end

% Update the new image and traces, and set the current slice & phase:
updateSlice(handles,sj,pj)




% ------------------------------------------------------------------------
function [hf,scrollAmount] = parseInputs(h,event)
% Switch for Matlab or Java Callbacks, or on manual call event==[]:
if isempty(event)
    hf = h;
    scrollAmount = 0;
    
elseif isa(event,'java.awt.event.MouseWheelEvent')
    % Using Java
    scrollAmount = event.getWheelRotation(); %event.getUnitsToScroll

    % Find the Matlab figure handle:
    set(0,'ShowHiddenHandles','on');
    kids = get(0,'Children');
    hf = findobj(kids,'-depth',0,'Name',h.getSpecifiedTitle.toCharArray');
    
else
    % Using Matlab
    scrollAmount = event.VerticalScrollCount;
    hf = h;
end
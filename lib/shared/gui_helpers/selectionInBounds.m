function selectionInBounds(~,eventdata)
%SELECTIONINBOUNDS Listener callback which keeps a list selection within
%   the bounds of the list items.
%
% Attach to an object with the following call:
%
%   addlistener(hOjbect,'String','PostSet',@selectionInBounds);
%
%

hobj = get(eventdata,'AffectedObject'); % Item we're listening to
selected = get(hobj,'Value');
contents = get(eventdata,'NewValue');

% Check current selection value
if isempty(selected)
    sel = 1;
else
    sel = selected(end);
end

% Keep selection within limits of the list items:
if isempty(contents)      
    sel = 1;                   % Empty default == 1
elseif sel > numel(contents)   
    sel = numel(contents);
elseif sel < 1
    sel = 1;
    %fprintf(2,'Shouldn''t have gotten here, what went wrong?\n');
    %keyboard
end

% Update the uicontrol, if necessary
if isempty(selected) || ~isequal(sel, selected)
    set(hobj,'Value',sel)
end


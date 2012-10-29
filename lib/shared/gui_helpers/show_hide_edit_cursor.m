function show_hide_edit_cursor(hObject)
%SHOW_HIDE_EDIT_CURSOR Show or hide the edit cursor depending on the
% checked state of hObject
%
% HOBJECT is the handle to the currently executing callback, which should
% be a uimenu item used to toggle the state of the edit cursor.
%
% This funcion is mainly used for development/debugging.

hf = ancestor(hObject,'figure');                % Figure
hcsr = findall(hf,'tag','Standard.EditPlot');   % Cursor uitoggletool
switch get(hcsr,'Visible')
    case 'on'
        set(hcsr,'Visible','off')
        set(hObject,'Checked','off')
    case 'off'
        set(hcsr,'Visible','on')
        set(hObject,'Checked','on')
end
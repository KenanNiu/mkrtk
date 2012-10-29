function toggleCheck(hObject)
%TOGGLECHECK Toggle the check on UIMENU item.
assert(ishandle(hObject))
assert(isequal(get(hObject,'Type'),'uimenu'),'TOGGLECHECK only works on UIMENU objects')

switch get(hObject,'Checked')
    case 'off'
        set(hObject,'Checked','on')
    case 'on'
        set(hObject,'Checked','off')
end
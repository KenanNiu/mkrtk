function tf = ischecked(hObject)
%ISCHECKED Check to see if UIMENU item is checked or not
assert(ishandle(hObject))
assert(isequal(get(hObject,'Type'),'uimenu'),'ISCHECKED only works on UIMENU objects')

switch get(hObject,'Checked')
    case 'on'
        tf = true;
    case 'off'
        tf = false;
end
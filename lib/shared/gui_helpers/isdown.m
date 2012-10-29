function tf = isdown(hObject)
%ISDOWN Check to see if UITOGGLETOOL is down or not
assert(ishandle(hObject))
assert(isprop(hObject,'State'),...
    'ISDOWN only works on UITOGGLETOOL objects, or objects with ''State'' property')

switch get(hObject,'State')
    case 'on'
        tf = true;
    case 'off'
        tf = false;
end


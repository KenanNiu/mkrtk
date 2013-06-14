function curs = getpointer(fig)
%GETPOINTER Complimentary function to SETPOINTER
if isappdata(fig,'Pointer')
    curs = getappdata(fig,'Pointer');
else
    curs = get(fig,'Pointer');
end
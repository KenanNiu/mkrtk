function krf2d(hfig,~)

% Key release function simply sets the 'CurrentKey' of the figure's appdata
% to an empty string
setappdata(hfig,'CurrentKey','');

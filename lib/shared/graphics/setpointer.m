function setpointer(fig,curs,cdata,hotspot)
%SETPOINTER Set the figure pointer
%
% We have this function for two reasons, namely
%   1. Mac doen't natively support the 'fleur' pointer
%   2. To record the name of custom pointers provided by setptr
%
% Usage:
%   setpointer(fig,cursor_name)             % cursor_name can be anything handled
%                                           % natively, or provided by setptr
%
%   setpointer(fig,'custom',cdata)          % set custom cursor
%   setpointer(fig,'custom',cdata,hotspot)  % set custom cursor with hotspot
%
% Examples:
%   setpointer(fig,'fleur')
%
%   setpointer(fig,'eraser')
%   getpointer(fig)
%       ans = 
%           'eraser'
%
%   setpointer(fig,'custom',rand(16),[8 8])

% In this function, the follow definitions apply:
%   POINTER - name that we will store in appdata
%   CURS    - name that will be passed to set() or setptr()

if ismac && strcmpi(curs,'fleur') % special case
    pointer = curs;         
    curs = 'custom';        
    cdata = fleurcursor();
    hotspot = [8 8];
else
    
    if strcmpi(curs,'custom')
        % setpointer(fig,'custom',cdata,[hotspot])
        pointer = curs;
        if ~exist('hotspot','var')
            hotspot = get(fig,'PointerShapeHotSpot');
        end
    else
        % setpointer(fig, cursor_name)
        std = ['crosshair', 'fullcrosshair', 'arrow', 'ibeam', 'watch', 'topl',...
            'topr', 'botl', 'botr', 'left', 'top', 'right', 'bottom', 'circle',...
            'cross', 'fleur', 'custom', 'hand'];
    
        if ~any( strcmpi(curs,std) )
            % curs is a setptr pre-defined custom name
            pointer = curs;
        end
    
    end
end

% Finally set the cursor
if strcmpi(curs,'custom')
    % manual set
    set(fig,'Pointer',curs,'PointerShapeCData',cdata,'PointerShapeHotSpot',hotspot)
else
    % use setptr for standard and standard special names:
    setptr(fig,curs)
end

% Record the pointer name:
setappdata(fig,'Pointer',pointer)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
function C = fleurcursor

C = [...
 NaN NaN NaN NaN NaN NaN   2   1   2 NaN NaN NaN NaN NaN NaN NaN
 NaN NaN NaN NaN NaN   2   1   1   1   2 NaN NaN NaN NaN NaN NaN
 NaN NaN NaN NaN   2   1   1   1   1   1   2 NaN NaN NaN NaN NaN
 NaN NaN NaN   2   2   2   2   1   2   2   2   2 NaN NaN NaN NaN
 NaN NaN   2   2 NaN NaN   2   1   2 NaN NaN   2   2 NaN NaN NaN
 NaN   2   1   2 NaN NaN   2   1   2 NaN NaN   2   1   2 NaN NaN
   2   1   1   2   2   2   2   1   2 NaN NaN   2   1   1   2 NaN
   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2
   2   1   1   2   2   2   2   1   2   2   2   2   1   1   2 NaN
 NaN   2   1   2 NaN NaN   2   1   2 NaN NaN   2   1   2 NaN NaN
 NaN NaN   2   2 NaN NaN   2   1   2 NaN NaN   2   2 NaN NaN NaN
 NaN NaN NaN   2   2   2   2   1   2   2   2   2 NaN NaN NaN NaN
 NaN NaN NaN NaN   2   1   1   1   1   1   2 NaN NaN NaN NaN NaN
 NaN NaN NaN NaN NaN   2   1   1   1   2 NaN NaN NaN NaN NaN NaN
 NaN NaN NaN NaN NaN NaN   2   1   2 NaN NaN NaN NaN NaN NaN NaN
 NaN NaN NaN NaN NaN NaN NaN   2 NaN NaN NaN NaN NaN NaN NaN NaN];

C( C==2 ) = NaN;


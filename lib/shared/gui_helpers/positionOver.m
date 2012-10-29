function positionOver(hOver,hUnder,opt)
%POSITIONGUI Position the gui hOver over the top of the calling gui hUnder
%
% POSITIONOVER(H,HUNDER) positions the figure H over the top of the figure
% HUNDER, centred laterally and slightly set in from the top.
%
% POSITIONOVER(H,HUNDER,'CENTER') positions the figure H centered laterally
% and vertically on top of HUNDER.

if nargin == 1 || ~isempty(hUnder)
    
    % Handle units
    cu = get(hUnder,'Units');
    fu = get(hOver,'Units');
    set([hUnder hOver],'Units','pixels')
    
    % Create correct position
    pc = get(hUnder,'Position');
    pf = get(hOver,'Position');
    x = pc(1) + (pc(3)-pf(3))/2;
    if nargin == 3 && isequal(opt,'center')
        y = pc(2) + (pc(4)-pf(4))/2;
    else
        y = pc(2) + pc(4) - 30 - pf(4);
    end
    pf(1:2) = [x y];    
    
    % Update position
    set(hOver,'Position',pf)
    
    % Restore units:
    set(hUnder,'Units',cu)
    set(hOver,'Units',fu);
    
else
    movegui(hOver,'center')
end
function [d,P1,P2] = momentArm(p1,p2,p3,p4)
%MOMENTARM Calculate moment arm from two lines, each defined by two points.
%
% See also distanceLines3d.m in geom3d toolbox

line1 = edge2line([p1 p2]);
line2 = edge2line([p3 p4]);

if nargout == 1
    % Marginally quicker if we don't require the points:
    d = mutual_perpendicular(line1,line2);
else
    [d,P1,P2] = mutual_perpendicular(line1,line2);
end

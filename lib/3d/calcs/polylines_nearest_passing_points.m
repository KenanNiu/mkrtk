function [p1,p2] = polylines_nearest_passing_points(poly1,poly2)
%POLYLINES_NEAREST_PASSING_POINTS Find the closest points that lie on each
% poly line.
% 
% Inputs:
%   POLY1 and POLY2 are each M-by-3 and N-by-3 lists of points in [x y z]
%   that describe two lines (polylines).
%
% Outputs:
%   P1 and P2 are 1-by-3 vectors of [x, y, z] describing the location of
%   the nearest point on each polyline.

% See also NEAREST_POINT_ON_POLYLINE

assert(size(poly1,2)==3 && size(poly2,2)==3, ...
    'Polylines must be N-by-3')

% Handle degenerate cases:
if (size(poly1,1) == 1) && (size(poly2,1) == 1)
    % Trivial case:
    p1 = poly1;
    p2 = poly2;
    return
    
elseif (size(poly1,1) == 1)
    p1 = poly1;
    [p2] = nearest_point_on_polyline(p1,poly2);
    return
    
elseif (size(poly2,1) == 1)
    p2 = poly2;
    [p1] = nearest_point_on_polyline(p2,poly1);
    return
    
end

% Create set of N-1 line segments, or edges: [x1 y1 z1, x2 y2 z2]
edges1 = [poly1(1:end-1,:) poly1(2:end,:)];
edges2 = [poly2(1:end-1,:) poly2(2:end,:)];

n1 = size(edges1,1);
n2 = size(edges2,1);

E1 = repmat(edges1,n2,1);
tmp = repmat(edges2',n1,1);
E2 = reshape(tmp(:),6,[])';

% Brute force, check all combinations:
[D,P1,P2] = distBetween2Segments(E1,E2);

% Find closest points and their segments:
[~,id] = min(D);
p1 = P1(id,:);
p2 = P2(id,:);


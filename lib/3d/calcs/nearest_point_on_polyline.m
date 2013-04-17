function [pt,di,idx,ti] = nearest_point_on_polyline(point, polyline)
%NEAREST_POINT_ON_POLYLINE Find the point 
% Inputs:
%   POINT is N-by-2 or N-by-3 (2D or 3D)
%   POLYLINE is N-by-2 or N-by-3 list of points which describe a curve [x y z]
%
% Outputs:
%   PT  is 1-by-2 or 1-by-3 point on the polyline which is closest to POINT
%   D   is the distance from POINT to PT 
%   IDX is the index of the segment of POLYLINE that PT lies on
%   T   is the parameterized length (0->1) along the segment, such that 
%           pt = p0 + ti*v
%        where p0 is the point:
%           p0 = polyline(idx,:) 
%        and v is the vector:
%           v = polyline(idx+1,:) - polyline(idx,:)

ndim = size(point,2);
assert(size(polyline,2) == ndim, ...
    'POINT & POLYLINE must have the same number of dimensions')
assert( (ndim == 2) || (ndim == 3), ...
    'Number of dimensions must be either 2 or 3')

% Convert polyline into list of segments:
seg = [polyline(1:end-1,:) polyline(2:end,:)];

% Get edge direction vector v, and parameterize the infinite line forms of
% seg in the form:  p0 + t*v
if ndim == 2
    v = seg(:,3:4) - seg(:,1:2);    % Edge direction vector:
    p0 = seg(:,1:2);
else
    v = seg(:,4:6) - seg(:,1:3);    % Edge direction vector:
    p0 = seg(:,1:3);
end

% Vector from line origin p0 to point
dp = bsxfun(@minus, point, p0);

% Parametrized distance from p0 along the line vector
t = bsxfun(@rdivide, sum(bsxfun(@times, dp, v), 2), sum(v.^2,2));

% Constrain projections to within bounds of segments (ie, 0 < t < 1)
t( t < 0 ) = 0;
t( t > 1 ) = 1;

% Coordinates of projected points:
if ndim == 2
    pp = bsxfun(@plus, p0, [t.*v(:,1) t.*v(:,2)]);
else
    pp = bsxfun(@plus, p0, [t.*v(:,1) t.*v(:,2) t.*v(:,3)]);
end

% Distance between point and its projection onto the segment:
if ndim == 2
    dist = sqrt( (point(:,1) - pp(:,1)).^2 + (point(:,2) - pp(:,2)).^2 );
else
    dist = sqrt( (point(:,1) - pp(:,1)).^2 + (point(:,2) - pp(:,2)).^2 + (point(:,3) - pp(:,3)).^2);
end

[di,idx] = min(dist);
pt = pp(idx,:);
ti = t(idx);
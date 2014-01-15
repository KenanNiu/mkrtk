function varargout = nearest_point_on_edge(point,edge)
% See DISTANCEPOINTEDGE3D.  This function operates much the same but with
% different outputs
%
% Usage:
% [pp,t,dist] = nearest_point_on_edge(point,edge)

% direction vector of each edge
vl = edge(:, 4:6) - edge(:, 1:3);

% compute position of points projected on the supporting line
% (Size of t is the max number of edges or points)
t = linePosition3d(point, [edge(:,1:3) vl]);

% ensure degenerated edges are correclty processed (consider the first
% vertex is the closest)
delta = vectorNorm3d(vl);
t(delta < eps) = 0;

% change position to ensure projected point is located on the edge
t(t < 0) = 0;
t(t > 1) = 1;

% coordinates of projected point
pp = bsxfun(@plus, edge(:,1:3), [t .* vl(:,1) t .* vl(:,2) t .* vl(:,3)]);

% process output arguments
varargout{1} = pp;
if nargout > 1
    varargout{2} = t;
    
    % compute distance between point and its projection on the edge
    dist = sqrt( sum( (point-pp).^2 ,2) );
    
    varargout{3} = dist;
end
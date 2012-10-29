function varargout = plotcloud(xyz,varargin)
% Plot a N-by-3 point cloud with plot options
h = plot3(xyz(:,1),xyz(:,2),xyz(:,3),varargin{:});
if nargout > 0
    varargout{1} = h;
end

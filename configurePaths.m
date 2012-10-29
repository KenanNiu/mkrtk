function configurePaths(mode)
% Set up the paths for the specified
%
% MODE is an optional string that can be one of the following
%   'all'       Set all paths
%   '2d'        Set paths for 2d Segmentation program
%   '3d'        Set paths for 3d Registration program

if nargin == 0 || isempty(mode)
    mode = 'all';
end

% Helper function
add = @(p)addpath(genpath(p));

% Generic paths:
add('lib/shared')

% Add 2D paths
if any( strcmpi(mode, {'all','2d'}) )
    add('lib/2d');
end

% Add 3D paths
if any( strcmpi(mode, {'all','3d'}) )
    add('lib/3d');
end
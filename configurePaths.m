function configurePaths(modestr)
% Set up the paths for the specified program.
%
% CONFIGUREPATHS when called by itself will add all necessary paths for all
% programs.  
%
% CONFIGUREPATHS(MODE) adds the paths for the specified run mode, which can
% be one of the following:
%   'all'       Set all paths
%   '2d'        Set paths for 2d Segmentation program
%   '3d'        Set paths for 3d Registration program

if nargin == 0 || isempty(modestr)
    modestr = 'all';
end

% Helper function
add = @(p)addpath(genpath(p));

% Generic paths:
add('lib/shared')

% Add 2D paths
if any( strcmpi(modestr, {'all','2d'}) )
    add('lib/2d');
end

% Add 3D paths
if any( strcmpi(modestr, {'all','3d'}) )
    add('lib/3d');
end

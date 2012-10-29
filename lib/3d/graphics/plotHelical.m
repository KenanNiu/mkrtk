function varargout = plotHelical(hax,axs,pt,varargin)
%PLOTHELICAL Plot helical axis.
%
% PLOTHELICAL(H,AXS,PT) plots the helical axis described by the 1-by-3
% direction vector AXS and the 1-by-3 xyz point PT onto axes H.
%
% PLOTHELICAL(H,AXS,PT,SCALE) plots the direction vector with length SCALE.
%
% PLOTHELICAL(...,Parameter1,Value1,Parameter2,Value2,...) passes the
% parameter/value pairs to the plotting functions.  Note that the
% parameters have to be compatible with QUIVER3 and PLOT3.
%
% HA = PLOTHELICAL(...) returns the handle to the axis' hggroup.
%
%   See also HELICALAXIS, SCREWAXIS.

% Joshua Martin, 02-Feb-2012

% ----- Defaults:
scale = 40;
linewidth = 1.6;
colour = 'b';


% ----- Parse specific inputs:

% Scale:
if ~isempty(varargin) && isnumeric(varargin{1})
    scale = varargin{1};
    varargin(1) = [];
else
    
end

% Line width:
lw_idx = find(strcmpi(varargin,'linewidth'));
if ~isempty(lw_idx)
    linewidth = varargin{lw_idx+1};
    varargin(lw_idx:lw_idx+1) = [];
end
 
% Colour:
c_idx = find(strcmpi(varargin,'color'));
if ~isempty(c_idx)
    colour = varargin{c_idx+1};
    varargin(c_idx:c_idx+1) = [];
end

% All other parameters get passed directly.



% ----- Plot the hggroup:

oldnp = get(hax,'NextPlot');
set(hax,'NextPlot','add');      % Set "hold on"

hsa = hggroup('Parent',hax,'HitTest','off');

% Bacground arrow (black)
quiver3(pt(1),pt(2),pt(3),axs(1),axs(2),axs(3),scale,...
    'parent',hsa,...
    'Color','k',...
    'LineWidth',linewidth,...
    varargin{:})
% Foreground arrow (colour)
quiver3(pt(1),pt(2),pt(3),axs(1),axs(2),axs(3),scale,...
    'parent',hsa,...
    'Color',colour,...
    'LineWidth',0.7*linewidth,...
    'LineStyle','-.',...
    varargin{:})
% Origin:
plot3(pt(1),pt(2),pt(3),'o',...
    'Parent',hsa,...
    'LineWidth',0.5*linewidth,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',colour,...
    'MarkerSize',linewidth*5,...
    varargin{:});

set(hax,'NextPlot',oldnp)   % Restore hold state

% ----- Assign output(s): 
if nargout > 0
    varargout{1} = hsa;
end

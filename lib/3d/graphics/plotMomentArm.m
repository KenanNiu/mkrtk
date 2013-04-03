function varargout = plotMomentArm(hax,pt1,pt2,varargin)
%PLOTMOMENTARM Plot moment arm defined by two points
%
%
%   See also MOMENTARM, CALCMOMENTARMS, PLOTHELICAL

% ----- Defaults:
linewidth = 1.6;
colour = 'r';

% ----- Parse specific inputs:

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

% Plot the hggroup:

hgg = hggroup('Parent',hax,'HitTest','off');

%oldnp = get(hax,'NextPlot');
%set(hax,'NextPlot','add');      % Set "hold on"

% Data:
X = [pt1(1) pt2(1)];
Y = [pt1(2) pt2(2)];
Z = [pt1(3) pt2(3)];

% Background line (black)
line(X,Y,Z,...
    'parent',hgg,...
    'Color','k',...
    'Linewidth',linewidth,...
    varargin{:});

% Foreground line (colour)
line(X,Y,Z,...
    'parent',hgg,...
    'Color',colour,...
    'Linewidth',linewidth,...
    'Linestyle','-.',...
    varargin{:});

% Endpoint circles
plot3(X,Y,Z,'o',...
    'Parent',hgg,...
    'LineWidth',0.5*linewidth,...
    'MarkerEdgeColor',colour,...
    'MarkerSize',linewidth*5,...
    varargin{:});

% Endpoint crosses
plot3(X,Y,Z,'x',...
    'Parent',hgg,...
    'Color','k',...
    'LineWidth',0.7*linewidth,...
    'MarkerSize',linewidth*5,...
    varargin{:});


%set(hax,'NextPlot',oldnp)   % Restore hold state

% ----- Assign output(s): 
if nargout > 0
    varargout{1} = hgg;
end

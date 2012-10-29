function varargout = LineSpec(varargin)
%LINESPEC Return line specifications for standard Matlab line 
%
% MAP = LINESPEC() returns a 5-by-2 cell array containing standard matlab
% LineSytle names in the first column and their corresponding line
% specifiers in the second column.
%
% LNAMES = LINESPEC(LSPEC) takes in a cell array of line specifiers and
% returns a cell array of LineStyle names.
%
% LSPEC = LINESPEC(LNAMES) takes in a cell array of LineStyle names and
% returns a cell array of line specifiers.
%
% Only minimal error checking.
%
%   See also COLOURSPEC.

% Joshua Martin, 10-Feb-2012

% Line definitions
map = {...            
    'solid','-';
    'dash','--';
    'dot',':';
    'dash-dot','-.';
    'none','none'};

if nargin == 0
    varargout{1} = map;
    return
end

arg1 = varargin{1};

% Determine the form of the input:
tf = false(size(map));
for j = 1:numel(arg1)
    tf = tf + strcmpi(map,arg1{j});
end

[~,c] = max(sum(tf));

% C now tells us which column in MAP the input corresponds to.

% So we define the output column as the opposite:
if c == 1
    k = 2;
elseif c == 2
    k = 1;
end

% Create output list and populate with default:
list = cell(size(arg1));        % empty list
[list{:}] = deal(map{1,k});     % default

% Now fill with mapped values
for j = 1:numel(arg1)
    tf = strcmpi(map(:,c),arg1{j});
    if any(tf)
        list{j} = map{tf,k};
    end
end

varargout{1} = list;

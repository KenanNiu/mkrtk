function varargout = current(prop,h,varargin)
%CURRENTPROP Central place for setting & retrieving current navigation
% states:
%   - Current Slice
%   - Current Phase
%
% Set usage:
%   current('slice',axeshandle,j)
%   current('phase',axeshandle,j)
%
% Get usage:
%   k = current('slice',axeshandle)
%   k = current('phase',axeshandle)
%
%
% Supercedes CURRENTMODE, CURRENTSLICE, CURRENTPHASE

% 16-Mar-2012, Joshua Martin


% Check the user has provided valid property & correct handle:
switch lower(prop)
    case 'slice'
        assert(isequal(get(h,'Type'),'axes'),'Must be handle to the axis system')
        
    case 'phase'
        assert(isequal(get(h,'Type'),'axes'),'Must be handle to the axis system')
        
    otherwise
        error('Unknown property %s',prop)    
end

% Create the property tag, one of:
%   - CurrentSlice
%   - CurrentPhase
tag = ['Current' upper(prop(1)) prop(2:end)];

% Run the action:
switch numel(varargin)
    case 0  % GET
        varargout{1} = getappdata(h,tag);
        
    case 1  % SET
        setappdata(h,tag,varargin{1});
        varargout = {};
end

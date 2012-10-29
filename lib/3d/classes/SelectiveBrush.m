classdef SelectiveBrush < handle
% SelectiveBrush object
%
%
% b = SelectiveBrush(h)
% b = SelectiveBrush([h1,h2,..])
%
% b.enable
% b.disable
% b.isenabled
%


properties (SetAccess = protected)
    axis
    targets
    brushtool
end % properties

properties (Hidden = true, SetAccess = protected)
    context_menu
    original_cbk
    children
    child_hvis
end

methods
    function b = SelectiveBrush(targHandles)   % Constructor
        if nargin == 1
            assert( all( ishghandle(targHandles) ), ...
                message('SelectiveBrush:SelectiveBrush:BadInputs'))
            
            assert( all( strcmpi(get(targHandles,'type'),'line') ), ...
                message('SelectiveBrush:SelectiveBrush:BadInputs'))
            
            p = ancestor(targHandles,'axes');
            assert( allequal(p), ...
                message('SelectiveBrush:SelectiveBrush:ParentDiffers')) 
            
            ax  = ancestor(targHandles(1),'axes');
            fig = ancestor(ax,'figure');
            b.axis = ax;
            b.targets = ( targHandles(:) )';
            b.brushtool = get_brushtool_handle(fig);
            b.original_cbk = get(b.brushtool,'ClickedCallback');
        else
            error(message('SelectiveBrush:SelectiveBrush:IncorrectNargin'))
        end
        
    end %SelectiveBrush()
    
    % ====================================================================
    
    % ------------------------------------
    function clearAllBrushing(b)
        hitem = findobj(b.context_menu,'Tag','BrushSeriesContextMenuClearAllBrushing');
        cbk = get(hitem,'Callback');
        % Version compatability:
        if iscell(cbk)                  % 2010b seems to have a cell array but 
            [fun, ax] = deal(cbk{:});   % function handle & axes handle
            feval( fun, hitem, [], ax ) % These are the inputs it requires
        else
            feval( cbk, hitem );        % Later versions just have the function handle
        end                 
    end %clearAllBrushing()
    
    % ------------------------------------
    function disable(b)
        if ~b.isenabled
            disp('Not enabled')
            return
        end
        b.clearAllBrushing();                                   % Remove all brushing
        set(b.brushtool,'State','off')                          % Toggle tool off
        feval(get(b.brushtool,'ClickedCallback'))               % Evaluate programatic callback
        set(b.brushtool,'ClickedCallback',b.original_cbk)       % Restore interactive callback
        try
            % This will fail if the children have changed (usually if the
            % plot has been re-drawn)
            set(b.children,{'HandleVis'},b.child_hvis)
        catch ME
            disp(ME)
        end
        b.children = [];
    end %disable()
    
    % ------------------------------------
    function enable(b)
        if b.isenabled
            disp('Already enabled')
            return
        end
        % First, turn handle visibility off for non-target objects
        %   (Axis handle visibiliyt needs to remain on though)
        b.children = findall(b.axis,'-not','Type','axes');   % findall would return b.axis in the result
        b.child_hvis = get( b.children, 'Visible');
        set(b.children(~ismember(b.children,b.targets)),'HandleVis','off')
        
        % Now configure and fire the databrushing tool
        cbk = @()putdowntext('brush',b.brushtool);  % Create new callback for programatic access
        set(b.brushtool,'ClickedCallback',cbk)      % Set callback
        set(b.brushtool,'State','on')               % Push the button
        feval(get(b.brushtool,'ClickedCallback'))   % Fire its callback
        
        % Store handle for context menu, which doesn't exist until
        % databrshing is turned on: 
        fig = ancestor(b.axis,'figure');
        b.context_menu = findall(fig,'tag','BrushSeriesContextMenu');
        configureActions(b.context_menu);
        
    end %enable()
    
    % ------------------------------------
    function D = getData(b)
        %GETDATA Get return the data for the brushed object(s)
        % Output is a cell array.
        n = numel(b.targets);
        D = cell(1,n);
        for j = 1:n
            ht = b.targets(j);
            x = get(ht,'XData');
            y = get(ht,'YData');
            z = get(ht,'ZData');
            D{j} = [x(:) y(:) z(:)];
        end
    end %getData()
        
    % ------------------------------------
    function tf = isenabled(b)
        tf = ~isequal(get(b.brushtool,'ClickedCallback'),b.original_cbk);
    end %isenabled()
    
    % ------------------------------------
    function setColor(b,clr)
        %SETCOLOR Set brushing colour
        fig = ancestor(b.axis,'figure');
        brushmode = getuimode(fig,'Exploration.Brushing');
        msd = get(brushmode,'ModeStateData');
        msd.color = clr;
        set(brushmode,'ModeStateData',msd)
    end %setColor()
    
end %methods
    

% Static methods
methods (Static)
    
    % ------------------------------------
    function HideTool(fig)
        %HIDETOOL Hide the brushing tool in the figure toolbar
        db = get_brushtool_handle(fig);
        set(db,'Visible','off')
    end %HideTool()
    
    % ------------------------------------
    function ShowTool(fig)
        %SHOWTOOL Show the brushing tool in the figure toolbar
        db = get_brushtool_handle(fig);
        set(db,'Visible','on')
    end %ShowTool
    
end %methods (Static)


end %classdef SelectiveBrush

% ------------------------------------------------------------------------
function tf = allequal(x)
%ALLEQUAL Test if all elements of x are equal
if numel(x) == 1
    tf = true;
    return
elseif iscell(x)
    tf = all( cellfun(@isequal,x,repmat(x(1),size(x))));
else
    tf = all(x == x(1));
end
end %allequal()


% ------------------------------------------------------------------------
function configureActions(context_menu)
items = get(context_menu,'children');
% Remove the "Replace with" action:
delete( findobj(items,'Tag','BrushSeriesContextMenuReplaceWith') )
end %configureActions()

% ------------------------------------------------------------------------
function h = get_brushtool_handle(fig)
% GET_BRUSHTOOL_HANDLE Get the handle to the brush tool on the specified
% figure
h = findall( fig, 'tag', 'Exploration.Brushing');
end

% ------------------------------------------------------------------------
function varargout = message(msgid,varargin)
%MESSAGE Generate warning or error messaged from message id
switch msgid
            
    case 'SelectiveBrush:SelectiveBrush:IncorrectNargin'
        msgstr = 'Incorrect number of input arguments';
    case 'SelectiveBrush:SelectiveBrush:BadInputs'
        msgstr = 'Inputs must be valid line object graphics handles'; 
    case 'SelectiveBrush:SelectiveBrush:ParentDiffers'
        msgstr = 'Graphics objects must have the same parent axes';
    
        
    otherwise
        msgstr = sprintf('Unhandled Message identifier: %s',msgid);
end

switch nargout
    case 1
        % Error message:
        msg.message = msgstr;
        msg.identifier = msgid;
        varargout = {msg};
        
    case 2
        % Warning message:
        varargout = {msgid,msgstr};
end

end %message()
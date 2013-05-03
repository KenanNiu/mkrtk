function phaseDisplay(handles)
% Display clouds in current phase from ICP results
%
% 

% Shorthand:
ha = handles.axes1;
hf = handles.figure1;

% Clear axes
delete(get(ha,'Children'))

% Return if nothing loaded:
if isempty(handles.Models)
    return
end

% Tool states:
STATIC  = isequal('on', get(handles.ShowStaticTool, 'State'));
DYNAMIC = isequal('on', get(handles.ShowDynamicTool,'State'));
CS      = isequal('on', get(handles.ShowAxesTool,   'State'));

% Retrieve some variables: 
p = get(handles.PhaseSlider,'Value');   % Current phase
nm = numel(handles.Models);             % Number of objects
clrs = get(ha,'ColorOrder');

% Motion:
[q,x] = models2qx(handles.Models,quaternion(NaN),NaN);


% Version specific plot opts:
%   7.12.0 --> 2011a
%   7.13.0 --> 2011b
cprops = {'Marker','.','LineStyle','none'};            % Default cloud plotting properties
qprops = {};                        % Quiver plotting properties
persistent LATE
if isempty(LATE)
    if verLessThan('matlab','7.12.0')  % 2011a and later
        LATE = false;
    else
        % 2011a and later
        LATE = true;
    end
end
if LATE
    cprops = [cprops,{'MarkerSize',8}];
    qprops = {'LineWidth',2.6};
end

     


% ---------- Draw models
for mj = 1:nm
    
    % Static models
    if STATIC && ~isempty(handles.Models(mj).HiRes) && handles.Models(mj).HiRes.Visible
        if ~isempty(q) && q(p,mj).isnan
            % If it's a NaN, this indicates that motion has been
            % calculated, but not for this position.  So it's appropriate
            % that the models are not draw for this case
        else
            if isempty(q)
                % In this case, no motion has been calculated at all, so
                % just display the static cloud in its initial position:
                cld = handles.Models(mj).HiRes;
            else
                % Motion has been calculated, so apply it and show the
                % model in its registered location:
                qj = q(p,mj);               % get quaternion
                xj = squeeze(x(p,mj,:))';   % orient displacement to 1-by-3
                cld = handles.Models(mj).HiRes.transform( qj, xj );
            end
            hc = cld.plot(ha,'Color',clrs(mj,:), cprops{:});
            set(hc,'uicontextmenu',handles.hscmenu) % Add menu 
            if CS
                cld.plotcs(ha)
            end
        end
        
    end
    
    % Dynamic models
    if DYNAMIC && ~isempty(handles.Models(mj).LoRes) && handles.Models(mj).LoRes(1).Visible
        cld = handles.Models(mj).LoRes(p);
        cld.plot(ha,'Color', [1 1 1]*0.1,'Linestyle','-',cprops{:});
        if CS
            cld.plotcs(ha)
        end
    end
    
end

% ---------- Draw Helical Axes
if ~isempty(handles.HelicalAxis)
    for j = 1:numel(handles.HelicalAxis)
        if handles.HelicalAxis(j).Visible && ~isempty(handles.HelicalAxis(j).Axis)
            
            c2id = strcmp({handles.Models.Tag},handles.HelicalAxis(j).Item2);
            hclr = clrs(c2id,:);
            
            % Plot single Helical Axis:
            plotHelical(ha,...
                handles.HelicalAxis(j).Axis(p,:),...
                handles.HelicalAxis(j).Point(p,:),...
                'Color','m',qprops{:})
            
            % Plot alternate Helical Axis
%             plotHelical(ha,...
%                 handles.HelicalAxis(j).Axis(p-1,:),...
%                 handles.HelicalAxis(j).Point(p-1,:),...
%                 'Color',hclr,qprops{:})

            % Plot all Helical Axes
%             for k = 4:size(handles.HelicalAxis(j).Axis,1)
%                 try
%                     plotHelical(ha,...
%                         handles.HelicalAxis(j).Axis(k,:),...
%                         handles.HelicalAxis(j).Point(k,:),...
%                         'Color','k',qprops{:})
%                 catch ME
%                     fprintf('Axis %2d of %2d not found...\n',k,size(handles.HelicalAxis(j).Axis,1))
%                 end
%             end

        end
    end
end

% ---------- Draw Line of Action
if ~isempty(handles.LineOfAction)
    for j = 1:numel(handles.LineOfAction)
        lj = handles.LineOfAction(j);
        if lj.Visible && isfield(lj,'Point') && isfield(lj,'Vector')
            % Draw it as the same as the helical axis
            plotHelical(ha, lj.Vector(p,:), lj.Point(p,:), 'Color', 'c', qprops{:})
        end
    end
end

% ---------- Draw Moment Arm
if ~isempty(handles.MomentArm)
    for j = 1:numel(handles.MomentArm)
        mj = handles.MomentArm(j);
        if mj.Visible && isfield(mj,'Point1') && isfield(mj,'Point2')
            plotMomentArm(ha,mj.Point1(p,:),mj.Point2(p,:),qprops{:});
        end
    end
end


% ----------- Configure slider & annotation
% Slider & annotation are visible if any one (or more) of the following are
% true:
%   - dynamic models are being displayed
%   - motion has been calculated (so static models move)
%nphases = max( [ cellfun(@numel,{handles.Models(:).LoRes}) 0] );
HASMOTION = any( ~cellfun(@isempty,{handles.Models(:).q}) );
an3dVis = 'off';
swf = [];
kpf = [];
if ( DYNAMIC || HASMOTION )    % dynamic models are displayed OR motion exists
    an3dVis = 'on';
    swf = @swf3d;
    kpf = @kpf3d;
end
set([handles.PhaseSlider,handles.PhaseText],'Visible',an3dVis)

% Now straddle the update to the swf & kpf with turning off of the
% exploration tools, then reverting them
explorationStates(handles.figure1,'off')    % Turn the tools off to avoid conflict
set(hf,'WindowScrollWheelFcn',swf)          %\_ Update functions
set(hf,'WindowKeyPressFcn',kpf)             %/
explorationStates(handles.figure1,'revert') % Revert the tools to user's current states

axis(ha,'tight')    % Tighten
zoom(ha,'reset')    % Set these axis limits as default


% ------------------------------------------------------------------------
function explorationStates(hf,opt)
%EXPLORATIONSTATES 
%
% There is a conflict between Matlab's internal WindowKeyPressFcn (and
% others?) and the functions that we use in this program when any of the
% exploration tools are down and we try to set our custom
% WindowKeyPressFcn.
%
% This function has two calls which should be used together to (first) set
% the exploration buttons to State='off', then the custom WindowKeyPressFcn
% can be added to the figure, then revert the exploration button states to
% their former values. An example would be:
%
%   explorationStates(hfig,'off')
%   set(hfig,'WindowKeyPressFcn',@mykpf)
%   explorationStates(hfig,'revert')
% 

persistent h

switch opt
    case {'off','up'}
        h = [0 0 0];
        if isdown(findall(hf,'Tag','Exploration.Rotate'))
            rotate3d(hf,'off');
            h(1) = 1;
        end
        if isdown(findall(hf,'Tag','Exploration.Pan'))
            pan(hf,'off');
            h(2) = 1;
        end
        if isdown(findall(hf,'Tag','Exploration.ZoomIn')) || ...
                isdown(findall(hf,'Tag','Exploration.ZoomOut'))
            set(0,'CurrentFigure',hf)   % Ensure hf is the current figure
            zoom % toggle zoom
            h(3) = 1;
        end
        
    case 'revert'
        if (h(1))
            rotate3d(hf,'on');
        end
        if (h(2))
            pan(hf,'on');
        end
        if (h(3))
            set(0,'CurrentFigure',hf)   % Ensure hf is the current figure
            zoom % toggle zoom
        end
                
        
    otherwise
        error('unknown option: %s',opt)
end



function draggable(hobj,varargin)
%
%
%   draggable(h)
%   draggable(h,'constrainto',hgobject)
%
%
% NOTE:
%   This function depends on the following:
%       - NEAREST_POINT_ON_POLYLINE
%       - MUTUAL_PERPENDICULAR
%       - GEOM3D toolbox (FEX)
%       - and probably others...

% Handle vector inputs:
if numel(hobj) > 1
    for j = hobj(:)'
        draggable(j,varargin{:})
    end
    return
end

[object_data] = parse_inputs(varargin{:});

fig = ancestor(hobj,'figure');
ax  = ancestor(hobj,'axes');

% Store all the figure data that we over-ride:
initial.wbmf = get(fig,'WindowButtonMotionFcn');
initial.wbuf = get(fig,'WindowButtonUpFcn');
initial.wbdf = get(fig,'WindowButtonDownFcn');
setappdata(fig,'draggable_initial_data',initial);

%set(hobj,'ButtonDownFcn',@button_down)

% Configure the object's control point
configure(hobj,object_data);

% Add properties
setappdata(hobj,'draggable_object_data',object_data)

% Enable draggable by setting the ButtonDownFcn on the object:
set(hobj,'ButtonDownFcn',@button_down)

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function hgt = configure(obj,opts)

ax = ancestor(obj,'axes');

% Configure the hgtransform object:
hgt = hgtransform('Parent',ax);
set(obj,'Parent',hgt)
setappdata(hgt,'target_object',obj)
setappdata(obj,'transform_object',hgt)

constraint      = opts.constraint_type;
constraint_data = opts.constraint_data;

% Object data:
xo = get(obj,'XData'); xo = xo(:);
yo = get(obj,'YData'); yo = yo(:);
zo = get(obj,'ZData'); zo = zo(:);

switch lower(constraint)
    
    case 'none'
        % Set the control point
        cp = mean( [xo, yo, zo], 1);
        manipulate(hgt,'setcontrolpoint',cp)
        
    case 'hghandle'
        hc = constraint_data;
        xc = get(hc,'XData');
        yc = get(hc,'YData');
        zc = get(hc,'ZData');
        switch get(hc,'type')
            
            case 'line'
                assert( isequal( get(obj,'type'), 'line'), 'hgobject type not supported')
                
                % OBJ is a line, and we're constraining it to HC which is
                % also a line.
                
                % Find nearest passing points:                
                [pc,po] = polylines_nearest_passing_points([xc(:) yc(:) zc(:)],[xo yo zo]);
                
                % Set the control point on the object:
                manipulate(hgt,'setcontrolpoint',po)
                
                % Move object to the point on the constraining object:
                manipulate(hgt,'moveto',pc)                
                
                
            otherwise
                error('unhandled hg type')
        end
        
    otherwise
        error('unhandled constraint type')
        
end %switch


% ------------------------------------------------------------------------
function button_down(obj,~)
fig = gcf;
ax  = gca;
set(fig,'WindowButtonMotionFcn',{@motion_fcn,obj})
set(fig,'WindowButtonUpFcn',{@button_up,obj})
hgt = getappdata(obj,'transform_object');
set_stored_current_point(hgt,get(ax,'CurrentPoint'));    % Initialize from axes


% ------------------------------------------------------------------------
function button_up(fig,~,obj)

% Re-set properties that we have temporarily overridden:
initial = getappdata(fig,'draggable_initial_data');
set(fig,'WindowButtonMotionFcn',initial.wbmf)
set(fig,'WindowButtonUpFcn',initial.wbuf)
%set(fig,'WindowButtonDownFcn',initial.wbdf)

objdata = getappdata(obj,'draggable_object_data');

% Evaluate buttonupfcn:
if isa(objdata.btnupfcn,'function_handle')
    feval(objdata.btnupfcn,obj)
end


% ------------------------------------------------------------------------
function motion_fcn(fig,~,obj)
ax = ancestor(obj,'axes');

% Current cursor location
current_point_new = get(ax,'CurrentPoint');

% Graphics transform object:
hgt = getappdata(obj,'transform_object');
opts = getappdata(obj,'draggable_object_data');
constraint = opts.constraint_type;
constraint_data = opts.constraint_data;

switch lower(constraint)
    
    case 'none'
        % Free movement:
        
        % Last recorded cursor location
        current_point_old  = get_stored_current_point(hgt);    % Current axes point, 2-by-3
        
        % Find displacement:
        dp = mean( current_point_new - current_point_old, 1 );
        
        % Translate object:
        manipulate(hgt,'translate',dp)
        
        
    case 'hghandle'
        % In this case a graphics object is specified for the constraint
        hc = constraint_data;
        
        % Handle each differty type of constraint:
        switch get(hc,'Type')
            case 'line'
                curve_xyz = [get(hc,'XData')',get(hc,'YData')',get(hc,'ZData')'];
                
                control_point_new = constrained_displacement_line(current_point_new,curve_xyz);
                
                manipulate(hgt,'moveto',control_point_new)

                                
            otherwise
                error('unhandled graphics object type')
        end
    
        %case {'linesegment','line_segment'}
        % Constrained by list of 2D or 3D points which define a line
                
            
    otherwise
        error('unhandled constraint')
        
end

% Update the stored current point
set_stored_current_point(hgt,current_point_new)


% ------------------------------------------------------------------------
function cp_3d_coinc = constrained_displacement_line(current_point,curve_xyz)
% SELECTED_POINT should be in the format of the 'CurrentPoint' axes
% property
%
%
p0 = current_point(1,:);
nml = normalizeVector3d( diff(current_point) );

plane = createPlane(p0,nml);

crv_proj = projPointOnPlane(curve_xyz,plane); % Project onto plane
crv_2d = planePosition(crv_proj,plane);     % Get planar (u,v) coordinates
cp_2d  = planePosition(p0,plane);           % Just in case

[~,~,cp_seg_id,cp_t] = nearest_point_on_polyline(cp_2d, crv_2d);

% Now convert back into 3d location using curve parameters:
pt_on_curve = @(xyz,seg,t) xyz(seg,:) + t*(xyz(seg+1,:) - xyz(seg,:));
cp_3d_coinc = pt_on_curve(curve_xyz,cp_seg_id,cp_t);%curve_xyz(cp_seg_id,:) + cp_t*(curve_xyz(cp_seg_id+1,:) - curve_xyz(cp_seg_id,:));



% ------------------------------------------------------------------------
function [objdata] = parse_inputs(varargin)
% opt.constraint_type       [ {'none'} | 'hghandle' ]
% opt.constraint_data       [  { [] }  | 'hghandle' ]
p = inputParser;
addParamValue(p,'constrainto',[]);
addParamValue(p,'buttonupfcn',[]);
addParamValue(p,'endfcn',[])
parse(p,varargin{:})

% Configure constraint
objdata.constraint_type = 'none';
objdata.constraint_data = p.Results.constrainto;
if ~isempty(objdata.constraint_data)
    if ishghandle(objdata.constraint_data)
        objdata.constraint_type = 'hghandle';
    end
end

% Record buttonupfcn:
objdata.btnupfcn = p.Results.buttonupfcn;

% Configure endfcn
objdata.endfcn = p.Results.endfcn;



% ------------------------------------------------------------------------
function pt = get_stored_current_point(obj)
pt = getappdata(obj,'AxesCurrentPoint');

% ------------------------------------------------------------------------
function set_stored_current_point(obj,pt)
assert( numel(pt) == 6 )
setappdata(obj,'AxesCurrentPoint',pt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Constraint functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------------------
function [pointc] = constrain_with_line_segment(point,data)

if ishghandle(data)
    h = data;
    lx = get(h,'XData');
    ly = get(h,'YData');
    lz = get(h,'ZData');
    line_xyz = [lx(:),ly(:),lz(:)];
else
    line_xyz = data;
end

% Insure we have full 3D specification:
if size(line_xyz,2) < 3
    line_xyz(:,3) = 0;
end

pointc = nearest_point_on_polyline(point, line_xyz);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Geometry functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



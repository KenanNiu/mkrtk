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

% Store all the data that we over-ride:
initial.wbmf = get(fig,'WindowButtonMotionFcn');
initial.wbuf = get(fig,'WindowButtonUpFcn');
initial.wbdf = get(fig,'WindowButtonDownFcn');
setappdata(fig,'draggable_initial_data',initial);

% Configure the object's control point
configure_cp(hobj,object_data)

% Add properties
%setappdata(fig,'draggable_shared_data',shared_data)
setappdata(hobj,'draggable_object_data',object_data)

% Enable draggable by setting the ButtonDownFcn on the object:
set(hobj,'ButtonDownFcn',@button_down)


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function configure_cp(obj,opts)

constraint      = opts.constraint_type;
constraint_data = opts.constraint_data;

switch lower(constraint)
    
    case 'none'
        % do nothing
        
    case 'hghandle'
        hc = constraint_data;
        switch get(hc,'type')
            
            case 'line'
                assert( isequal( get(obj,'type'), 'line'), 'hgobject type not supported')
                
                % OBJ is a line, and we're constraining it to HC which is
                % also a line.
                % >> Set control point as closest point:
                xc = get(hc,'XData');
                yc = get(hc,'YData');
                zc = get(hc,'ZData');
                xo = get(obj,'XData');
                yo = get(obj,'YData');
                zo = get(obj,'ZData');
                [pc,po] = polylines_nearest_passing_points([xc(:) yc(:) zc(:)],[xo(:) yo(:) zo(:)]);
                
                d3 = pc - po;
                xo = xo + d3(1);
                yo = yo + d3(2);
                zo = zo + d3(3);
                
                set(obj,'XData',xo)
                set(obj,'YData',yo)
                set(obj,'ZData',zo)
                
                set_control_point(obj,pc)
                
                % In this case
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
set_stored_current_point(obj,get(ax,'CurrentPoint'));    % Initialize from axes


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

current_point_old  = get_stored_current_point(obj);    % Current axes point, 2-by-3
current_point_new = get(ax,'CurrentPoint');

opts = getappdata(obj,'draggable_object_data');
constraint = opts.constraint_type;
constraint_data = opts.constraint_data;

% Get point data of the object:
x = get(obj,'XData');
y = get(obj,'YData');
z = get(obj,'ZData');


switch lower(constraint)
    
    case 'none'
        % Free movement:
        % Does not use the control point (draggable_ControlPoint)
        
        dp = mean( current_point_new - current_point_old, 1 );    % Displacement
        x = x + dp(1);      %\
        y = y + dp(2);      % > New location
        z = z + dp(3);      %/
        
        set_stored_current_point(obj,current_point_new);
        
    case 'hghandle'
        % In this case a graphics object is specified for the constraint
        hc = constraint_data;
        
        % These objects also use a control point to handle movement:
        control_point_old = get_control_point(obj);
        
        switch get(hc,'Type')
            case 'line'
                curve_xyz = [get(hc,'XData')',get(hc,'YData')',get(hc,'ZData')'];
                
                control_point_new = constrained_displacement_line(control_point_old,current_point_new,curve_xyz);
                
                % Set new control point:
                dp = control_point_new - control_point_old;
                
                % New location:
                x = x + dp(1);
                y = y + dp(2);
                z = z + dp(3);
                                
            otherwise
                error('unhandled graphics object type')
        end
        set_control_point(obj,control_point_new);
    
        %case {'linesegment','line_segment'}
        % Constrained by list of 2D or 3D points which define a line
                
            
    otherwise
        error('unhandled constraint')
        
end

set(obj,'XData',x,'YData',y,'ZData',z)




% ------------------------------------------------------------------------
function cp_3d_coinc = constrained_displacement_line(previous_point,current_point,curve_xyz)
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


return
% Now to find the displacement, we need to do the same thing for 
% previous_point
pp_2d = planePosition(previous_point(1,:),plane);
[~,~,pp_seg_id,pp_t] = nearest_point_on_polyline(pp_2d, crv_2d);
pp_3d_coinc = pt_on_curve(curve_xyz,pp_seg_id,pp_t);

dp = cp_3d_coinc - pp_3d_coinc;

%plot3(cp_3d_coinc(1),cp_3d_coinc(2),cp_3d_coinc(3),'mx')
%plot3(pp_3d_coinc(1),pp_3d_coinc(2),pp_3d_coinc(3),'gx')


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
function pt = get_control_point(obj)
pt = getappdata(obj,'draggable_ControlPoint');

% ------------------------------------------------------------------------
function set_control_point(obj,pt)
assert( numel(pt) == 3 )
setappdata(obj,'draggable_ControlPoint',pt)

% ------------------------------------------------------------------------
function pt = get_stored_current_point(obj)
pt = getappdata(obj,'draggable_AxesCurrentPoint');

% ------------------------------------------------------------------------
function set_stored_current_point(obj,pt)
assert( numel(pt) == 6 )
setappdata(obj,'draggable_AxesCurrentPoint',pt);

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



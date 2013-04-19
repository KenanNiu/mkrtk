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
hgt = configure(hobj,object_data);

% Add properties
setappdata(hgt,'draggable_object_data',object_data)

% Enable draggable by setting the ButtonDownFcn on the object:
set(hobj,'ButtonDownFcn',@button_down)

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function hgt = configure(obj,opts)

ax = ancestor(obj,'axes');

% Configure the hgtransform object:
hgt = hgtransform('Parent',ax);
set(obj,'Parent',hgt)

% Store handles:
setappdata(hgt,'target_object',obj)
setappdata(obj,'transform_object',hgt)

% Get options:
constraint      = opts.constraint_type;
constraint_data = opts.constraint_data;
allow_rotate    = opts.allowrotate;

% Object data:
xo = get(obj,'XData'); xo = xo(:);
yo = get(obj,'YData'); yo = yo(:);
zo = get(obj,'ZData'); zo = zo(:);

plot3d = @(xyz,varargin)plot3(xyz(:,1),xyz(:,2),xyz(:,3),varargin{:});

if allow_rotate
    
    % Build locations for control handles:
    xyzo = [xo yo zo];      
    p0 = mean( xyzo, 1 );   
    
    if numel(xo) == 4
        % Use points on edges
        cpoints = diff( [xyzo; xo(1) yo(1) zo(1)], 1)./2 + xyzo;
    else
        % Use vertices
        cpoints = xyzo;
    end
    
    % Add control handles as children of hgtransform
    
    props = {'Parent',hgt,'LineStyle','none','Marker','o','Color','k','MarkerFaceColor','c','MarkerSize',11};
    
    hrcp = plot3d(cpoints,props{:},'Tag','_RotateControl');    % 'Rotate' control points
    htcp = plot3d(p0,props{:},'Tag','_TranslateControl');      % 'Translate' control points
    
    % Enable draggable on the controls:
    set([hrcp,htcp],'ButtonDownFcn',@button_down)
        
end

% Now configure the object's control point:
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
% Here OBJ could be the target object, or any one of its controls, all of
% which are children of the hgtransform object.
fig = gcf;
ax  = gca;
hgt = ancestor(obj,'hgtransform');  % Should be immediate parent, but just in case

cp = get(ax,'CurrentPoint');        % Initialize from axes
set_stored_current_point(hgt,cp);   % Store

% In some cases we need to know which data point in OBJ was selected on
% buttondown, so find it:
v = normalizeVector3d( diff( cp, 1 ) );              % View vector
x = get(obj,'XData'); x = x(:);
y = get(obj,'YData'); y = y(:);
z = get(obj,'ZData'); z = z(:);  
d = distancePointLine3d( [x y z], [cp(1,:) v] );
[~,id] = min(d);

% Store the id:
setappdata(hgt,'SelectedObject',obj);
setappdata(hgt,'SelectedPointIndex',id);
    
% Finally configure the motion function & button up function:
set(fig,'WindowButtonMotionFcn',{@motion_fcn,obj})
set(fig,'WindowButtonUpFcn',{@button_up,obj})


% ------------------------------------------------------------------------
function button_up(fig,~,obj)

% Re-set properties that we have temporarily overridden:
initial = getappdata(fig,'draggable_initial_data');
set(fig,'WindowButtonMotionFcn',initial.wbmf)
set(fig,'WindowButtonUpFcn',initial.wbuf)

hgt = ancestor(obj,'hgtransform');  % Should be immediate parent, but just in case
objdata = getappdata(hgt,'draggable_object_data');

% Evaluate user's buttonupfcn:
if isa(objdata.btnupfcn,'function_handle')
    feval(objdata.btnupfcn,obj)
end


% ------------------------------------------------------------------------
function motion_fcn(fig,~,obj)
% OBJ is a handle to the graphics object we're manipulating.  It could
% either be the primary object of interest, or one of the controls that are
% built for manipulation.

hgt = ancestor(obj,'hgtransform');
ax = ancestor(hgt,'axes');

% Current cursor location
current_point_new = get(ax,'CurrentPoint');

% Object data & options:
opts = getappdata(hgt,'draggable_object_data');
constraint = opts.constraint_type;
constraint_data = opts.constraint_data;
allow_rotate = opts.allowrotate;


% -------------- Handle rotation commands ----------------- %
if allow_rotate && strcmpi(get(obj,'Tag'),'_RotateControl')
    
    % Last recorded cursor location
    current_point_old  = get_stored_current_point(hgt);    % Current axes point, 2-by-3
    
    % View vector:
    v = normalizeVector3d( diff( current_point_new, 1 ) );
    
    % Position of central control:
    hctr = findobj(hgt,'tag','_TranslateControl');
    p0 = [ get(hctr,'XData') get(hctr,'YData') get(hctr,'ZData') ];
    
    
    % Position of current control:
    rctrlxyz = [ get(obj,'XData')'  get(obj,'YData')'  get(obj,'ZData')' ];
    id = getappdata(hgt,'SelectedPointIndex');
    p1 = rctrlxyz(id,:);
    
    
    % We create a mean plane for orientation:
    C = princomp(rctrlxyz);
    mp_u = C(:,1)';
    mp_v = C(:,2)';
    nml = cross(mp_u,mp_v);
    
    % Vector to central control
    u = normalizeVector3d( p1 - p0 );
        
    % Rotation axis is the cross product of the mean plane normal with the
    % vector back to the central control:
    axs = cross(nml,u);
        
    % Calculate the angle of rotation
    pr = createPlane(p0,axs);   % in the plane of the rotation
    
    pt_old = intersectLinePlane( [current_point_old(1,:) v], pr );
    pt_new = intersectLinePlane( [current_point_new(1,:) v], pr );
    
    %ang = anglePoints3d(pt_old,p0,pt_new)
    dp1 = normalizeVector3d(pt_old-p0);
    dp2 = normalizeVector3d(pt_new-p0);
    ang = acos( dot(dp1,dp2,2) );
    
    if ~isreal(ang) % dp1 == dp2
        return
    end
    
    sign = @(x)(x>=0)*2-1;  % same as Matlab's sign(x), but assigns sign(0)==+1;
    %s = sign( dot(cross(nml, dp1), dp2 ) )
    n = normalizeVector3d( cross(dp1,dp2) );
    s = sign( dot( n, axs ) );
    %keyboard
    
    manipulate(hgt,'rotate','axisrotate',axs,s*ang)
    
    set_stored_current_point(hgt,current_point_new)
    return
    
end


% ------------ The code below is for translation ----------- %
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
% opt.allowrotate           [ {false} | true ]
p = inputParser;
addParamValue(p,'constrainto',[]);
addParamValue(p,'allowrotate',false);
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

% Record other parameters:
objdata.allowrotate = (p.Results.allowrotate == true) || any( strcmpi(p.Results.allowrotate,{'on','true'}) );
objdata.btnupfcn    = p.Results.buttonupfcn;
objdata.endfcn      = p.Results.endfcn;


% ------------------------------------------------------------------------
function pt = get_stored_current_point(hgt)
pt = getappdata(hgt,'AxesCurrentPoint');

% ------------------------------------------------------------------------
function set_stored_current_point(hgt,pt)
assert( numel(pt) == 6 )
setappdata(hgt,'AxesCurrentPoint',pt);

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



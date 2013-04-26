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
%
%   TODO:   - Mouse cursors - rotate/pan for control points
%           - Refactor
%           - Handle out-of-bounds rotations
%           - cache functions
%           - It seems that using the axis rotation tool drops all
%             draggable functionality.... why?
%           - Create a shunt for different calling options


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
initial.kpf  = get(fig,'KeyPressFcn');
initial.krf  = get(fig,'KeyReleaseFcn');

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

% Set the interaction mode
setappdata(fig,'draggable_InteractionMode','translate')

% And also configure the keypress functions for changing mode:
% rotate/translate
set(fig,'KeyPressFcn',@key_down)


% ------------------------------------------------------------------------
function button_up(fig,~,obj)

% Re-set properties that we have temporarily overridden:
initial = getappdata(fig,'draggable_initial_data');
set(fig,'WindowButtonMotionFcn',initial.wbmf)
set(fig,'WindowButtonUpFcn',initial.wbuf)

hgt = ancestor(obj,'hgtransform');  % Should be immediate parent, but just in case
objdata = getappdata(hgt,'draggable_object_data');

% Remove the rotation/translation mode control callback (only KeyPressFcn -
% KeyReleaseFcn is reset/removed at the end of its execution)
set(fig,'KeyPressFcn',initial.kpf)

% Now run The KeyReleaseFcn if it is still set:
krf = get(fig,'KeyReleaseFcn');
if isequal(krf,@key_up)
    % Create eventdata (key) structure:
    k = struct('Key',get(gcf,'CurrentKey'),'Character',get(gcf,'CurrentCharacter'));
    % Run callback:
    krf(fig,k)
end

% Evaluate user's buttonupfcn:
if isa(objdata.btnupfcn,'function_handle')
    feval(objdata.btnupfcn,obj)
end

% ------------------------------------------------------------------------
function R = crosshair_rotation_matrix(ax,current_point,p0)

[~,r] = spaceball(ax);

% The mouse cursor intersects the sphere centered at p0 here:
vhat = normalizeVector3d( diff( current_point, 1 ) );
vline = [current_point(1,:) vhat];      % view line
lsi = intersectLineSphere(vline, [p0 r]); %
ps = lsi(1,:);       % Should be the first point, if not, use linePosition3d to find front point

% Manipulation crosshairs are described by the intersection between
% the sphere and the two orthogonal planes that run from the the point
% ps through the sphere centre, p0
v_s0 = normalizeVector3d( p0 - ps );    % vector back to centre

R = vec_vec_rotationmatrix([1 0 0],v_s0);


% ------------------------------------------------------------------------
function motion_fcn(fig,~,obj)
% OBJ is a handle to the graphics object we're manipulating.  It could
% either be the primary object of interest, or one of the controls that are
% built for manipulation.

hgt = ancestor(obj,'hgtransform');
ax = gca; %ancestor(hgt,'axes');

interaction_mode = getappdata(fig,'draggable_InteractionMode');

% Current cursor location
current_point_new = get(ax,'CurrentPoint');

% Object data & options:
opts = getappdata(hgt,'draggable_object_data');
constraint = opts.constraint_type;
constraint_data = opts.constraint_data;
allow_rotate = opts.allowrotate;



if allow_rotate && isequal(interaction_mode,'rotate')
    
    % -------------- ROTATION -------------- %
    % 1. Rotate the object
    % 2. Rotate the LCS indicators
    % 3. Rotate the manipulation cross-hairs
    
    p0 = manipulate(hgt,'getcontrolpoint');
    
    current_point_old  = get_stored_current_point(hgt);
    
    R_old = crosshair_rotation_matrix(ax,current_point_old,p0);
    R_new = crosshair_rotation_matrix(ax,current_point_new,p0);
    
    dR = R_new/R_old;
    
    manipulate(hgt,'rotate',dR)
    
    hgt_lcs = findobj(ax,'Tag','rotation_LCS_destroy');
    manipulate(hgt_lcs,'rotate',dR)
    
    hgt_mxh = findobj(ax,'Tag','manipulation_crosshair_destroy');
    manipulate(hgt_mxh,'rotate',dR)
    
    
else
    
    
    % -------------- TRANSLATION -------------- %
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
    
end
% Update the stored current point
set_stored_current_point(hgt,current_point_new)


% ------------------------------------------------------------------------
function update_manipulation_crosshairs(ax,R,T)

% hgtransform tag:
TAG = 'manipulation_crosshair_destroy';

% Build transformation matrix:
H = [R T(:); 0 0 0 1];

% Check if hgtransform exists:
hgt_mxh = findobj(ax,'Type','hgtransform','Tag',TAG);

% Draw from scratch:
if isempty( hgt_mxh )
    
    
    [~,r] = spaceball(ax);
    
    theta = linspace(3*pi/2-pi/10,3*pi/2+pi/10,30);
    
    
    x1 = r*sin(theta);
    y1 = r*cos(theta);
    z1 = zeros(size(x1));
    
    x2 = r*sin(theta);
    y2 = zeros(size(x2));
    z2 = r*cos(theta);
    
    H = eye(4);
    H(1:3,1:3) = R;
    H(1:3,4) = T(:);
    
    % Configure the hgtransform:
    hgt_mxh = hgtransform('Parent',ax,'Matrix',H,'Tag',TAG);
    manipulate(hgt_mxh,'setcontrolpoint',T);
    
    % Plot the two curves:
    props = {'Color',[1 0.4 0],'Line','-','Parent',hgt_mxh};
    plot3(x1,y1,z1,props{:})
    plot3(x2,y2,z2,props{:})
    
else
    
    % Simple update:
    
    set(hgt_mxh,'Matrix',H)
end

% ------------------------------------------------------------------------
function update_LCS_indicators(ax,R,T)

% hgtransform tag:
TAG = 'rotation_LCS_destroy';

% Build transformation matrix:
H = [R T(:); 0 0 0 1];

% Check if hgtransform exists:
hgt_lcs = findobj(ax,'Type','hgtransform','Tag',TAG);

% Draw from scratch:
if isempty( hgt_lcs )
    
    xhat = [1 0 0];
    yhat = [0 1 0];
    zhat = [0 0 1];
    
    [~,r] = spaceball(ax);
    
    l = r/4;
    
    % vector2crosshair anonymous function:
    v2c = @(v) [-l*v; -l*v*.3; NaN NaN NaN; l*v*.3; l*v];
    
    % Configure the hgtransform:
    hgt_lcs = hgtransform('Parent',ax,'Matrix',H,'tag',TAG);
    manipulate(hgt_lcs,'setcontrolpoint',T)
    
    % Plot the three vectors:
    props = {'Color',[1 0.4 0],'Line','-','Parent',hgt_lcs};
    plot3d(ax,v2c(xhat),props{:});
    plot3d(ax,v2c(yhat),props{:});
    plot3d(ax,v2c(zhat),props{:});
    
else
    % Simple update:
    
    set(hgt_lcs,'Matrix',H)
    
end

% ------------------------------------------------------------------------
function key_down(fig,k)
%
% Flow:
%
%   Button Press (WindowButtonDownFcn)
%       '--> apply KeyPressFcn
%
%       Key Press
%           '--> apply KeyReleaseFcn
%
%       Key Release
%           '--> destroy KeyReleaseFcn
%
%       Button Release
%           '--> destroy KeyReleaseFcn
%
if strcmpi(k.Key,'shift') && isempty(k.Character)
    
    % ---------- Configure ---------- %
    set(fig,'KeyReleaseFcn',@key_up)
    
    disp('rotate mode')
    setappdata(fig,'draggable_InteractionMode','rotate')
    
    ax = gca;
    hgt = ancestor(gco,'hgtransform');  % Is there a better way to do this?
    p0 = manipulate(hgt,'getcontrolpoint');
    
    % ---------- Plot Spaceball bounds: ---------- %
    
    % Spaceball:
    [~,r] = spaceball( ax );
    
    % But in fact it is better to draw on the front of the bounding box, so
    % shift p0 toward the viewer along the view vector
    current_point = get(ax,'CurrentPoint');
    p0_front = p0 - diff( current_point, 1 )./2;
    
    % Create circle in (x,y) plane
    n = 100;
    tc = linspace(0,2*pi,n)';
    xc = sin(tc)*r;
    yc = cos(tc)*r;
    zc = zeros(size(xc));
    nml = [0 0 1];
    
    % Now transform into view plane
    xyz = rotate_to_view_plane([xc yc zc],nml,ax);
    xyz = xyz + p0_front(ones(n,1),:);
    
    % Plot:
    plot3d(ax,xyz,'Tag','_destroy','Color',[1 0.4 0],'Line','--')
    
    % ============= Note on the following: ============== %
    %   In what follows, the Matrix of the hgtransform will need to be
    %   updated at every mouse move
    
    hgt = ancestor(gco,'hgtransform'); % should probably do this a better way
    H = get(hgt,'Matrix');
    R = H(1:3,1:3);
    
    update_LCS_indicators(ax,R,p0_front)
    
    
    % ---------- Now plot manipulation crosshairs ---------- %
    
    R = crosshair_rotation_matrix(ax,current_point,p0);
    
    update_manipulation_crosshairs(ax,R,p0_front)
    
    
    
end


% ------------------------------------------------------------------------
function xyz = rotate_to_view_plane(xyzin, nml, ax)
% Transform the data XYZIN with normal NML so that NML lines up with the
% direction vector that the axes AX is being viewed from

% Get the two vectors:
v1 = nml;
current_point = get(ax,'CurrentPoint');
v2 = normalizeVector3d( diff(current_point) );   % View vector

% Now calcuate the rotation matrix which transforms v1 to v2:
R = vec_vec_rotationmatrix(v1,v2);

% Rotate
xyz = (R*xyzin')';


% ------------------------------------------------------------------------
function R = vec_vec_rotationmatrix(v1,v2)
% Calculate the rotation matrix that transforms v1 to v2, such that:
%   v2 = R*v1

% Get the axis and angle
angle = acos(v1*v2');
axs = cross(v1,v2)/norm(cross(v1,v2));

% A skew symmetric representation of the normalized axis
axis_skewed = [ 0 -axs(3) axs(2) ; axs(3) 0 -axs(1) ; -axs(2) axs(1) 0];

% Rodrigues formula for the rotation matrix
R = eye(3) + sin(angle)*axis_skewed + (1-cos(angle))*axis_skewed*axis_skewed;


% ------------------------------------------------------------------------
function [p0, r] = spaceball( ax )
% Create location & radius for space manipulation sphere

bounds = reshape(axis(ax), 2, []);

% Bounding circle, with centre
p0 = mean( bounds );
r = mean( diff(bounds./2) );    % Calculate a radius


% ------------------------------------------------------------------------
function key_up(fig,k)
% This funtion needs to be robust about destroying - items may or may not
% exist
if strcmpi(k.Key,'shift') && isempty(k.Character)
    setappdata(fig,'draggable_InteractionMode','translate')   % Revert to default []
    
    % Destroy all stuffs
    delete( findobj(fig,'-regexp','Tag','_destroy') )
end

% Automatically destroy / reset this callback:
initial = getappdata(fig,'draggable_initial_data');
set(fig,'KeyReleaseFcn',initial.krf)


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
function varargout = plot3d(xyz,varargin)
% Convenient wrapper function for plot3
%
% PLOT3D(XYZ)
% PLOT3D(XYZ,'PROPERTY',VALUE,...)
% PLOT3D(AX,...)

% Manage inputs:
if isscalar(xyz) && ishghandle(xyz)
    ax = xyz;
    xyz = varargin{1};
    varargin(1) = [];
    
    h = plot3(ax,xyz(:,1),xyz(:,2),xyz(:,3),varargin{:});
    
else
    
    h = plot3(xyz(:,1),xyz(:,2),xyz(:,3),varargin{:});
    
end

% Manage outputs
if nargout == 0
    varargout = {};
else
    varargout = {h};
end


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



function draggable(hobj,varargin)
%
%
%   draggable(h)
%   draggable(h,'ConstrainTo',hgobject)
%   draggable(h,'AllowRotate',true)
%
% NOTE:
%   This function depends on the following:
%       - NEAREST_POINT_ON_POLYLINE
%       - MUTUAL_PERPENDICULAR
%       - GEOM3D toolbox (FEX)
%       - and probably others...
%
%   TODO:   - Constraint stuff is a little under developed (or over-developed?)
%           - Refactor
%               - appdata
%           - Handle right-clicks during left-click for rotation (as an
%             alternative to shift)
%           - cache functions
%           - It seems that using the axis rotation tool drops all
%             draggable functionality.... why?
%           - Create a shunt for different calling options
%               - 'off' case
%
% 
% NOTE:
%   For cases using the 'ConstrainTo' option, it is a good idea to set the
%   'HitTest' of the constraint object to 'off' to prevent it from
%   capturing clicks.
%


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Object data:
xo = get(obj,'XData'); xo = xo(:);
yo = get(obj,'YData'); yo = yo(:);
zo = get(obj,'ZData'); zo = zo(:);


% Now configure the object's control point:
switch lower(constraint)
    
    case 'none'
        % Set the control point
        cp = mean( [xo, yo, zo], 1 );
        manipulate(hgt,'setcontrolpoint',cp)
        
    case 'hghandle'
        hc = constraint_data;
        xc = get(hc,'XData');
        yc = get(hc,'YData');
        zc = get(hc,'ZData');
        
        switch get(hc,'Type')
            
            case 'line'     % Constraint is a line object
                
                switch get(obj,'Type')
                    
                    case 'line'                 % Constraint==line, object==line
                        % Nearest points on object & constratin:
                        [pc,po] = polylines_nearest_passing_points([xc(:) yc(:) zc(:)],[xo yo zo]);
                        
                    case {'patch','surface'}    % Constraint==line, object==[ patch | surface ]
                        po = mean( [xo yo zo], 1 );                             % point on object
                        pc = nearest_point_on_polyline(po,[xc(:) yc(:) zc(:)]); % point on constraint 
                        
                    otherwise
                        error('Unhandled graphics object to constrain: %s',get(obj,'Type'))
                end
                % ---- Snap to line: ----
                
                % Set the control point on the object:
                manipulate(hgt,'setcontrolpoint',po)
                
                % Move object to the point on the constraining object:
                manipulate(hgt,'moveto',pc)
                
                % ----               ----
                
            otherwise
                error('Unhandled hg constraint type: %s',get(hc,'Type'))
            
        end %switch constraint hghandle type 
        
                
    otherwise
        error('Unhandled constraint type: %s',constraint)
        
end %switch constraint type


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Callback functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
function button_down(obj,~)
%BUTTON_DOWN Called when a mouse button is clicked on OBJ.
%
% The primary purpose of this function is to enable the
% WindowButtonMotionFcn, which does all the hard work of translating or
% rotating the object.
%

% Get handles:
fig = gcf;
ax  = gca;
hgt = ancestor(obj,'hgtransform');  % Should be immediate parent, but just in case

% Initialise the stored version of 'CurrentPoint', held in hgt's appdata:
cp = get(ax,'CurrentPoint');        % Initialize from axes
set_stored_current_point(hgt,cp);   % Store

% In some cases we need to know which data point in OBJ was selected on
% buttondown, so find it:
v = normalizeVector3d( diff( cp, 1 ) );             % View vector
x = get(obj,'XData'); x = x(:);                     %\
y = get(obj,'YData'); y = y(:);                     % > All the object's data
z = get(obj,'ZData'); z = z(:);                     %/

% Find the index of the point closest to the line defined by
% 'CurrentPoint':
d = distancePointLine3d( [x y z], [cp(1,:) v] );
[~,id] = min(d);

% Store the id:
setappdata(hgt,'SelectedObject',obj);
setappdata(hgt,'SelectedPointIndex',id);

% Set the interaction mode & appropriate mouse cursor
set_interaction_mode(fig,'translate')

% And also configure the keypress functions for changing mode (ie,
% rotate/translate):
set(fig,'KeyPressFcn',@key_down)

% Finally configure the motion function & button up function:
set(fig,'WindowButtonMotionFcn',{@motion_fcn,obj})
set(fig,'WindowButtonUpFcn',{@button_up,obj})


% ------------------------------------------------------------------------
function button_up(fig,~,obj)
%BUTTON_UP Called when a mouse button is release on OBJ
% 
% The primary purpose of this function is to disable the
% WindowButtonMotionFcn, and restore figure properties that were
% temporarily over-ridden for interaction

% Reset cursor:
setpointer(fig,'arrow');

% Re-set properties that we have temporarily overridden:
initial = getappdata(fig,'draggable_initial_data');
set(fig,'WindowButtonMotionFcn',initial.wbmf)
set(fig,'WindowButtonUpFcn',initial.wbuf)

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
hgt = ancestor(obj,'hgtransform');  % Should be immediate parent, but just in case
objdata = getappdata(hgt,'draggable_object_data');
if isa(objdata.btnupfcn,'function_handle')
    feval(objdata.btnupfcn,obj)
end


% ------------------------------------------------------------------------
function motion_fcn(fig,~,obj)
%MOTION_FCN Perform all the manipulation of OBJ on every mouse move event
%
% MOTION_FCN operates in interaction modes, 'rotate' or 'translate'.  This
% mode determines how OBJ is manipulated.

% Get the interaction mode:
interaction_mode = get_interaction_mode(fig);

% Get handles:
hgt = ancestor(obj,'hgtransform');
ax = gca; 

% Current cursor locations, old & new
current_point_old  = get_stored_current_point(hgt);    % Current axes point, 2-by-3
current_point_new = get(ax,'CurrentPoint');

% Object data & options:
opts = getappdata(hgt,'draggable_object_data');
constraint = opts.constraint_type;
constraint_data = opts.constraint_data;
allow_rotate = opts.allowrotate;


if allow_rotate && isequal(interaction_mode,'rotate')
    
    % -------------- ROTATION -------------- %
    
    %  In this mode we have to not only rotate the object, but also rotate
    %  the local coordinate indicators, as well as the manipulation
    %  cross-hairs
    
    % 1. Rotate the object
    % 2. Rotate the LCS indicators
    % 3. Rotate the manipulation cross-hairs
    
    % 0. Calculate the rotation required since last update:
    p0 = manipulate(hgt,'getcontrolpoint');
        
    r = spaceball(ax);
    
    [R_old,OOB_old] = crosshair_rotation_matrix(r,current_point_old,p0);
    [R_new,OOB_new] = crosshair_rotation_matrix(r,current_point_new,p0);
    
    dR = R_new/R_old;
    
    % Check for the in-plane rotation case when the cursor is outside the
    % bounds of the manipulation sphere;
    if OOB_old && OOB_new
        % Re-hash - we need to calculate an axial rotation around a vector
        % centred at p0 and in the direction of the current_point vectors.
        % We use only the component of the rotation in dR which occurs
        % about the axis aligned with current_point.
        
        % View axis:
        axs = normalizeVector3d( diff( current_point_new, 1 ) );
        
        [dRaxs,dRang] = rotmat2axisangle(dR);
        ang = dRang * dot(dRaxs,axs);
        H = makehgtform('axisrotate',axs,ang);
        dR = H(1:3,1:3);
        R_new = dR*R_old;
        
    else
        % Calculate simple change:
        dR = R_new/R_old;
        
    end
       
    
    % 1. Rotate the object    
    manipulate(hgt,'rotate',dR)
    
    % 2. Rotate the local coordinate system:
    p0_front = move_to_front(ax,p0);
    update_LCS_indicators(ax,R_new,p0_front)
    
    % 3. Rotate the manipulation cross-hairs:
    update_manipulation_crosshairs(ax,R_new,p0_front)
    
    % -------------- ROTATION -------------- %
    
else
    
    
    % -------------- TRANSLATION -------------- %
    switch lower(constraint)
        
        case 'none'
            % Free movement:
            
            
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
                    
                    %case {'linesegment','line_segment'}
                    % Constrained by list of 2D or 3D points which define a
                    % line, but are not plotted
                    
                otherwise
                    error('Unhandled graphics object type for constraint: %s', constraint)
                    
            end
            
            
            
            
        otherwise
            error('unhandled constraint')
            
    end
    % -------------- TRANSLATION -------------- %
    
    
end
% Update the stored current point
set_stored_current_point(hgt,current_point_new)


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

hgt = ancestor(gco,'hgtransform');  
opts = getappdata(hgt,'draggable_object_data');
allow_rotate = opts.allowrotate;

if allow_rotate && strcmpi(k.Key,'shift') && isempty(k.Character)
    
    % ---------- Configure ---------- %
    set(fig,'KeyReleaseFcn',@key_up)
    set_interaction_mode(fig,'rotate')
    
    ax = gca;
    
    p0 = manipulate(hgt,'getcontrolpoint');
    
    % ---------- Plot Spaceball bounds: ---------- %
    
    % Spaceball:
    r = spaceball( ax );
    
    % It is better to draw on the front of the bounding box, so
    % shift p0 toward the viewer along the view vector
    current_point = get(ax,'CurrentPoint');
    p0_front = move_to_front(ax,p0);
        
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
    
    R = crosshair_rotation_matrix(r,current_point,p0);
    
    update_manipulation_crosshairs(ax,R,p0_front)
    
    
end


% ------------------------------------------------------------------------
function key_up(fig,k)
% This function is run on keyboard key up event, or called by button_up
if strcmpi(k.Key,'shift') && isempty(k.Character)
    
    if numel(dbstack) == 1
        % Called by key_up event directly
        set_interaction_mode(fig,'translate')
    else
        % Called by button_up
        set_interaction_mode(fig,'off')
    end
    
    % Destroy all items tagged for destruction
    delete( findobj(fig,'-regexp','Tag','_destroy') )
    
    % Automatically destroy / reset this key_up callback:
    initial = getappdata(fig,'draggable_initial_data');
    set(fig,'KeyReleaseFcn',initial.krf)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
function [R,OOB] = crosshair_rotation_matrix(radius,current_point,p0)

% Cursor line: line made by current_point
A = current_point(1,:);
B = current_point(2,:);
vhat = normalizeVector3d( B-A );
cline = [A vhat];      % cursor line

% The cursor line should intersect the sphere; ie, d <= r
d = distancePointLine3d(p0,cline);

% Check out of bounds:
OOB = d > radius;
% If d > r, cline does not intersect the sphere.
if OOB    
    % Move cline so it is tangent to the sphere in order to ensure solution:
    tol = 1e-10;
    rhat = normalizeVector3d( (A + linePosition3d(p0,cline)*vhat) - p0 );
    cline(1:3) = p0 + ( rhat*(radius-tol) ); % -eps ensures we're never outside
end

% The sphere centered at p0 is pierced by the mouse cursor here:
lsi = intersectLineSphere(cline, [p0 radius]); %
ps = lsi(1,:);       % Should be the first point, if not, use linePosition3d to find front point

% Manipulation crosshairs are described by the intersection between
% the sphere and the two orthogonal planes that run from the the point
% ps through the sphere centre, p0
v_s0 = normalizeVector3d( p0 - ps );        % vector back to centre
R1 = vec_vec_rotationmatrix([1 0 0],vhat);  % crosshairs defined with x-basis; rotate to view line
R2 = vec_vec_rotationmatrix(vhat,v_s0);     % Rotate from view line to current line through sphere
R = R2*R1;



% ------------------------------------------------------------------------
function p_front = move_to_front(ax,point)
% For drawing manipulation controls in front of everything else.  perhaps
% not the best solution, an implementation is not formally correct, but
% it'll do for now... 
current_point = get(ax,'CurrentPoint');
p_front = point - diff( current_point, 1 )./2;


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
function [axs,ang] = rotmat2axisangle(R)
%ROTMAT2AXISANGLE Convert rotation matrix to axis / angle representation
%
% This function does the same job as vrrotmat2vec if you have simulink
% installed.  However, this function does not cater for the ambiguity when
% angle = pi, since we expect only small rotations to be passed in for this
% application.
tol = 1e-12;
mtc = trace(R);
if ( abs( mtc - 3 ) <= tol )
    axs = [1 0 0];
    ang = 0;
else
    ang = acos( (mtc - 1) / 2 );
    axs = [R(3,2)-R(2,3) R(1,3)-R(3,1) R(2,1)-R(1,2)] ./ (2*sin(ang)); 
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
    
    r = spaceball(ax);
    
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
function update_manipulation_crosshairs(ax,R,T)

% hgtransform tag:
TAG = 'manipulation_crosshair_destroy';

% Build transformation matrix:
H = [R T(:); 0 0 0 1];

% Check if hgtransform exists:
hgt_mxh = findobj(ax,'Type','hgtransform','Tag',TAG);

% Draw from scratch:
if isempty( hgt_mxh )
    
    
    r = spaceball(ax);
    
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
function [r] = spaceball( ax )
% Create location & radius for space manipulation sphere

bounds = reshape(axis(ax), 2, []);

% Bounding circle, with centre
%p0 = mean( bounds );
r = mean( diff(bounds./2) );    % Calculate a radius



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Utility functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
function curs = getpointer(fig)
%GETPOINTER Complimentary function to SETPOINTER
if isappdata(fig,'Pointer')
    curs = getappdata(fig,'Pointer');
else
    curs = get(fig,'Pointer');
end

% ------------------------------------------------------------------------
function setpointer(fig,curs,cdata,hotspot)
%SETPOINTER Set the figure pointer
%
% We have this function for two reasons, namely
%   1. Mac doen't natively support the 'fleur' pointer
%   2. To record the name of custom pointers provided by setptr
%
% Usage:
%   setpointer(fig,cursor_name)             % cursor_name can be anything handled
%                                           % natively, or provided by setptr
%
%   setpointer(fig,'custom',cdata)          % set custom cursor
%   setpointer(fig,'custom',cdata,hotspot)  % set custom cursor with hotspot
%
% Examples:
%   setpointer(fig,'fleur')
%
%   setpointer(fig,'eraser')
%   getpointer(fig)
%       ans = 
%           'eraser'
%
%   setpointer(fig,'custom',rand(16),[8 8])

% In this function, the follow definitions apply:
%   POINTER - name that we will store in appdata
%   CURS    - name that will be passed to set() or setptr()

if ismac && strcmpi(curs,'fleur') % special case
    pointer = curs;         
    curs = 'custom';        
    cdata = fleurcursor();
    hotspot = [8 8];
else
    
    if strcmpi(curs,'custom')
        % setpointer(fig,'custom',cdata,[hotspot])
        pointer = curs;
        if ~exist('hotspot','var')
            hotspot = get(fig,'PointerShapeHotSpot');
        end
    else
        % setpointer(fig, cursor_name)
        std = ['crosshair', 'fullcrosshair', 'arrow', 'ibeam', 'watch', 'topl',...
            'topr', 'botl', 'botr', 'left', 'top', 'right', 'bottom', 'circle',...
            'cross', 'fleur', 'custom', 'hand'];
    
        if ~any( strcmpi(curs,std) )
            % curs is a setptr pre-defined custom name
            pointer = curs;
        end
    
    end
end

% Finally set the cursor
if strcmpi(curs,'custom')
    % manual set
    set(fig,'Pointer',curs,'PointerShapeCData',cdata,'PointerShapeHotSpot',hotspot)
else
    % use setptr for standard and standard special names:
    setptr(fig,curs)
end

% Record the pointer name:
setappdata(fig,'Pointer',pointer)


% ------------------------------------------------------------------------
function pt = get_stored_current_point(hgt)
pt = getappdata(hgt,'AxesCurrentPoint');

% ------------------------------------------------------------------------
function set_stored_current_point(hgt,pt)
assert( numel(pt) == 6 )
setappdata(hgt,'AxesCurrentPoint',pt);


% ------------------------------------------------------------------------
function mode = get_interaction_mode(fig)
%GET_INTERACTION_MODE Complimentary function for SET_INTERACTION_MODE
TAG = 'draggable_InteractionMode';
mode = getappdata(fig,TAG);

% ------------------------------------------------------------------------
function set_interaction_mode(fig,mode)
%SET_INTERACTION_MODE For setting the interaction mode of the figure
TAG = 'draggable_InteractionMode';
switch mode
    case 'off'
        if isappdata(fig,TAG)
            rmappdata(fig,TAG)
        end
    case 'rotate'
        setappdata(fig,TAG,'rotate')
        setpointer(fig,'rotate')
    case 'translate'
        setappdata(fig,TAG,'translate')
        setpointer(fig,'fleur')
    otherwise
        error('unhandled interaction mode')
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Constraint functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
%                    Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
function C = fleurcursor

C = [...
 NaN NaN NaN NaN NaN NaN   2   1   2 NaN NaN NaN NaN NaN NaN NaN
 NaN NaN NaN NaN NaN   2   1   1   1   2 NaN NaN NaN NaN NaN NaN
 NaN NaN NaN NaN   2   1   1   1   1   1   2 NaN NaN NaN NaN NaN
 NaN NaN NaN   2   2   2   2   1   2   2   2   2 NaN NaN NaN NaN
 NaN NaN   2   2 NaN NaN   2   1   2 NaN NaN   2   2 NaN NaN NaN
 NaN   2   1   2 NaN NaN   2   1   2 NaN NaN   2   1   2 NaN NaN
   2   1   1   2   2   2   2   1   2 NaN NaN   2   1   1   2 NaN
   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2
   2   1   1   2   2   2   2   1   2   2   2   2   1   1   2 NaN
 NaN   2   1   2 NaN NaN   2   1   2 NaN NaN   2   1   2 NaN NaN
 NaN NaN   2   2 NaN NaN   2   1   2 NaN NaN   2   2 NaN NaN NaN
 NaN NaN NaN   2   2   2   2   1   2   2   2   2 NaN NaN NaN NaN
 NaN NaN NaN NaN   2   1   1   1   1   1   2 NaN NaN NaN NaN NaN
 NaN NaN NaN NaN NaN   2   1   1   1   2 NaN NaN NaN NaN NaN NaN
 NaN NaN NaN NaN NaN NaN   2   1   2 NaN NaN NaN NaN NaN NaN NaN
 NaN NaN NaN NaN NaN NaN NaN   2 NaN NaN NaN NaN NaN NaN NaN NaN];

C( C==2 ) = NaN;

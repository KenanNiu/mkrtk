function varargout = manipulate(obj,action,varargin)
% Graphics object manipulator

% Setter methods:
%   manipulate(hobj,'setcontrolpoint',pt)
%   manipulate(hobj,'translate',[dx dy dz])
%   manipulate(hobj,'translateto',[x y z])      
%   manipulate(hobj,'move',[dx dy dz])          % alias for 'translate'       
%   manipulate(hobj,'moveto',[x y z])           % alias for 'translateto'       
%   manipulate(hobj,'rotate',R)                 % rotate around control point with matrix R 
%   manipulate(hobj,'rotate',...)               % rotation with inputs as per hgtransform 
%
% Getter methods:
% pt = manipulate(hobj,'getcontrolpoint')


varargout = {};

% Run the requested function:
fname = lower(action);
fun = str2func( fname );

% Call appropriately according to nargout
n = nargout(fname);
if n == 0
    fun(obj, varargin{:});
else
    [varargout{1:n}] = fun(obj, varargin{:});
end


% ------------------------------------------------------------------------
function move(obj,dcp)

% Translate the object:
switch lower( get(obj,'type') )
    
    case 'hgtransform'
        T = get(obj,'Matrix');
        T(1:3,4) = T(1:3,4) + dcp(:);
        set(obj,'Matrix',T)
        
    otherwise
        x = get(obj,'XData') + dcp(1);
        y = get(obj,'YData') + dcp(2);
        z = get(obj,'ZData') + dcp(3);
        set(obj,'XData',x,'YData',y,'ZData',z)
end

% Adjust its control point:
cp_old = getcontrolpoint(obj);
cp_new = cp_old + dcp;
setcontrolpoint(obj,cp_new)


% ------------------------------------------------------------------------
function moveto(obj,cp_new)
cp_old = getcontrolpoint(obj);
move(obj, cp_new - cp_old);


% ------------------------------------------------------------------------
function rotate(obj,varargin) %#ok<DEFNU>
switch lower( get(obj,'type') )
    
    case 'hgtransform'
        H = get(obj,'Matrix');
        
        
        if isequal( size(varargin{1}), [3,3] )
            % manipulate(hobj,'rotate',R)
            R = varargin{1};
            Hr = [R, [0;0;0]; 0 0 0 1];
        else
            % manipulate(hobj,'rotate','xrotate',t)
            % manipulate(hobj,'rotate','axisrotate',axis,angle)
            % ... etc ...
            Hr = makehgtform(varargin{:});
        end
        cp = getcontrolpoint(obj);
        
        % Apply rotation to current Matrix
        H = Hr*H;
        
        % Now un-shift        
        cp2 = (Hr(1:3,1:3)*cp')';
        dt = cp - cp2;
        H(1:3,4) = H(1:3,4) + dt(:);
        
        % Update:
        set(obj,'Matrix',H)
                
        
    otherwise
        error('not implemented')
end

% Control point should not have changed


% ------------------------------------------------------------------------
function translate(obj,varargin) %#ok<DEFNU>
% Alias for 'move'
move(obj,varargin{:})

% ------------------------------------------------------------------------
function translateto(obj,varargin) %#ok<DEFNU>
% Alias for 'moveto'
moveto(obj,varargin{:})

% ------------------------------------------------------------------------
function cp = getcontrolpoint(obj)
cp = getappdata(obj,CPTAG);
assert( ~isempty(cp), ...
    'Control point has not been initialised. Use SETCONTROLPOINT to initialise.')

% ------------------------------------------------------------------------
function setcontrolpoint(obj,pt)
setappdata(obj,CPTAG,pt)

% ------------------------------------------------------------------------
function tag = CPTAG
tag = 'm_control_point';

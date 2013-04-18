function varargout = manipulate(obj,action,varargin)
% Graphics object manipulator

% manipulate(hobj,'initialize',[x y z])
% manipulate(hobj,'initialize',cobj)
% manipulate(hobj,'get')
% manipulate(hobj,'set',pt)
% manipulate(hobj,'move',[dx dy dz])


varargout = {};

% Run the requested function:
%feval( lower(action), obj, varargin{:} );
fname = lower(action);
fun = str2func( fname );

% Call appropriately according to nargout
n = nargout(fname);
if n == 0
    fun(obj, varargin{:});
else
    [varargout{1:n}] = fun(obj, varargin{:});
end


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


function moveto(obj,cp_new)
cp_old = getcontrolpoint(obj);
move(obj, cp_new - cp_old);


function translate(obj,varargin)
% Alias for 'move'
move(obj,varargin{:})

function translateto(obj,varargin)
% Alias for 'moveto'
moveto(obj,varargin{:})

function cp = getcontrolpoint(obj)
cp = getappdata(obj,CPTAG);
assert( ~isempty(cp), ...
    'Control point has not been initialised. Use SETCONTROLPOINT to initialise.')


function setcontrolpoint(obj,pt)
setappdata(obj,CPTAG,pt)


function tag = CPTAG
tag = 'm_control_point';

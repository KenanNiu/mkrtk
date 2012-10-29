classdef plane3d
%PLANE3D class
%
%
%


properties (SetAccess = protected)
    p0 = point3d();
    u  = vector3d();
    v  = vector3d();
end

properties (Dependent = true, SetAccess = private)
    n   % plane normal; also available through method: .normal
end

methods % ordinary methods
    
    function p = plane3d( varargin )
        clsIn = cellfun(@class,varargin,'UniformOutput',false);
        switch nargin
            case 0
                % default values
            case 1  % plane3d( point3d )
                pts = varargin{1};
                clsEx = {'point3d'};
                if isequal( clsIn, clsEx ) && ( numel(pts)==3 )
                    p = plane3d.throughpoints(pts(1),pts(2),pts(3));
                elseif ispoint
                    error('need exactly 3 points')
                elseif isnumeric(pts)
                    error('Creation from numveris not yet implemented')
                else
                    error('need 3x point3d')
                end
            case 2  % plane3d( point3d, vector3d )
                [pt, nml] = deal(varargin{:});
                clsEx = {'point3d', 'vector3d'};
                if isequal( clsIn, clsEx )
                    p = plane3d.pointnormal(pt,nml);
                elseif isnumeric(pt) && isnumeric(nml)
                    error('Creation from numerics not yet implemented')
                else
                    error(['Two inputs should be (poin3d, vector3d). ',...
                        'You provided (%s, %s).'],clsIn{:})
                end
            case 3  % plane3d( point3d, vector3d, vector3d )
                [pt,arg1,arg2] = deal(varargin{:});
                clspvv = {'point3d', 'vector3d', 'vector3d'};
                clsppp = {'point3d', 'point3d',  'point3d' };
                if isequal( clsIn, clspvv )
                    p = plane3d.pointdirections(pt,arg1,arg2);
                elseif isequal( clsIn, clsppp )
                    p = plane3d.throughpoints(pt,arg1,arg2);
                elseif isnumeric(pt) && isnumeric(arg1) && isnumeric(arg2)
                    error('Creation from numerics not yet implemented')
                else
                    error(['Three inputs should be (point3d, vector3d, vector3d). ',...
                        'You provided (%s, %s, %s).'],clsIn{:})
                end
                
            otherwise
                error('Incorrect number of input arguments')
        end %swtich nargin
    end %plane3d()
    
    % ----------------------------------------------------
    function tf = coincident(p,obj,tol)
        %COINCIDENT Test coincidence of two objects
        %
        % Plane-Plane:
        %   Planar coincidence/equality is assessed based on a 
        %   hierarchy two conditions:
        %       a) Plane normals are parallel (equal/equivalent)
        %       b) p2 origin lies in the plane of p1
        % Plane-Line:
        %   Coincidence is determined based on 
        %
        % Plane-Point:
        %   True if point lies on the plane
        if ~exist('tol','var')
            tol = eps;
        end
        siz = size(p);
        switch lower( class(obj) )
            case 'plane3d'
                p = p(:);       %\_ format convenience
                obj = obj(:);   %/
                tfp  = parallel(p,obj);             % Planes parallel
                tfc  = coincident(p,[obj.p0]',tol); % Plane2 origin coincident with Plane1
                tf = tfp & tfc;                     % Both must be true
                fprintf(2,'Can normals can be parallel but opposite?\n');
                keyboard
            case 'line3d'
                keyboard
            case 'point3d'
                d = p.distPointToPlane(obj);
                tf = chop(d) == 0;
            otherwise
                error('Unhandled object')
        end
        tf = reshape(tf,siz);
    end %coincident()
    
    % ----------------------------------------------------
    function display(p)
        %DISPLAY Called during unsupressed output to the 
        % command line.
        % See also disp() for standard un-formatted output.
        siz = size( p );
        nel = [1 cumprod( siz )];
        ndim = length(siz);
        disp(' ')
        fprintf('%s = \n\n',inputname(1));
        sstr = sprintf('%dx',siz);  sstr(end) = [];
        fprintf('  %s <a href="matlab:help plane3d">plane3d</a>',sstr);
        fprintf(' with <a href="matlab:methods(''plane3d'')">methods</a>\n\n');
        for iel = 1:nel(end)
            if nel(end) == 1
                sub = '';
            else
                sub = ')';
                jel = iel - 1;
                for idm = ndim : -1 : 1
                    idx = floor( jel / nel(idm) ) + 1;
                    sub = [',' int2str(idx) sub]; %#ok<AGROW>
                    jel = rem( jel, nel(idm) );
                end
                sub(1) = '(';
            end 
            pstr = sprintf('(%g, %g, %g)',p(iel).p0.xyz);
            ustr = sprintf('(%g, %g, %g)',p(iel).u.uvw);
            vstr = sprintf('(%g, %g, %g)',p(iel).v.uvw);
            fprintf( ['%s' sub ':\tp0: %s\tu: %s\tv: %s\n'],...
                inputname(1), pstr, ustr, vstr)
        end
        fprintf('\n')        
    end %display()
    
    % ----------------------------------------------------
    function d = distPointToPlane(p,pt)
        %DISTPOINTTOPLANE Normal distance of point to plane
        %   Tol is an optional 
        %
        % Uses Hessian form, ie : N.p = d
        % I this case, d can be found as : -N.p0, when N is normalized
        siz = size(p);
        normal_d = [p.normal.uvw]';   
        origin_d = [p.origin.xyz]';
        point_d  = [pt.xyz]';
        d = -sum(bsxfun(@times, normal_d, bsxfun(@minus, origin_d, point_d)), 2);
        d = reshape(d,siz);
    end %distPointToPlane()
    
    % ----------------------------------------------------
    function tf = eq( p1, p2 )
        %EQ Test equality of elements.
        % Usage:
        %   p1 == p2
        %   p1.eq(p2)
        % 
        % Testing equality is defined as testing coincidence:
        tf = coincident(p1,p2);        
    end %eq()
    
    % ----------------------------------------------------
    function nml = get.n( p )
        %N Get method to retrieve object normal
        nml = p.normal;
    end %n()
    
    % ----------------------------------------------------
    function i = intersect(p,obj)
        %INTERSECT Intersect object with plane
        % N planes and M objects not allowed.  Either:
        %   - N planes & N objects
        %   - N planes & 1 objects
        %   - 1 plane & N objects
        nump = numel(p);
        numo = numel(obj);
        if (nump ~= numo) && (min(nump,numo) > 1)
            error(['Inputs must have the same number of elements',...
                ' or one must be 1']);
        elseif nump > numo
            siz = size(p);
        else 
            siz = size(obj);
        end
        
        switch lower( class(obj) )
            case 'line3d'
                nml_d = [p.normal.uvw]';
                p_origin_d = [p.origin.xyz]';
                l_origin_d = [obj.origin.xyz]';
                l_vec_d    = [obj.direction.uvw]';
                % difference between origins of plane and line
                dp = bsxfun(@minus, p_origin_d, l_origin_d);
                % dot product of line direction with plane normal
                denom = sum(bsxfun(@times, nml_d, l_vec_d), 2);
                % relative position of intersection point on line (can be inf in case of a
                % line parallel to the plane)
                t = sum(bsxfun(@times, nml_d, dp),2) ./ denom;
                % compute coord of intersection point
                point = bsxfun(@plus, l_origin_d,  bsxfun(@times, [t t t], l_vec_d));
                i = point3d(point');
            otherwise
                error('Unhandled object')
        end
        i = reshape(i,siz);
        
    end %intersect()
    
    % ----------------------------------------------------
    function nml = normal( p )
        %NORMAL Get plane normal
        % Faster than using get.n
        uhat = [p.u];
        vhat = [p.v];
        nml = uhat.cross(vhat);
        nml = reshape( nml.normalize, size(p) );
    end
    
    % ----------------------------------------------------
    function tf = neq( p1, p2 )
        tf = ~eq(p1,p2);
    end %neq()
        
    % ----------------------------------------------------
    function p = normalize(p)
        u1 = [p.u];
        v1 = [p.v];
        u1 = u1.normalize;
        nml = cross(u1,v1);
        nml = nml.normalize;
        v1 = cross(nml,u1);       % v will now be normalized;
        for j = 1:numel(p)
            p(j).u = u1(j);
            p(j).v = v1(j);
        end
    end %normalize()
    
    % ----------------------------------------------------
    function o = origin(p)
        %ORIGIN Get plane origin
        o = reshape( [p.p0], size(p) );
    end %origin()
        
    % ----------------------------------------------------
    function tf = parallel( p1, p2 )
        %PARALLEL Test if planes are parallel
        % Parallel planes have equivalent normal vectors
        tf = [p1.n] == [p2.n];
        tf = reshape(tf,size(p1));
    end %parallel()
    
    % ----------------------------------------------------
    function varargout = plot(p,varargin)
        %PLOT Draw plane on current or specified axes
        % Usage:
        %   plane.plot
        %   PLOT(plane)
        %   PLOT(plane,axs)
        %   PLOT(plane,'PropertyName',propertyvalue,...)
        %   PLOT(plane,axs,'PropertyName',propertyvalue,...)
        %
        % Where properties are those accepted by PATCH
        if numel(p) ~= 1
            warning('Can we act on multiple planes? Loop?')
            keyboard
        end
        h = drawplane(p,varargin{:});        
        if nargout > 0
            varargout{1} = h;
        end
    end %plot()
end %methods

methods (Static)
    
    % ---------------------------------------------------
    function p = pointdirections(pt,uin,vin)
        %warn if norm of cross product is very small...
        % Create plane from point (point3d), and two orthogonal direction vectors
        % (vector3d)
        assert( all(numel(pt)== [numel(uin), numel(vin)]), ...
            'Number of elements must be the same')        
        p = plane3d.empty;
        for j = numel(pt) : -1 : 1
            p(j).p0 = pt(j);
            p(j).u  = uin(j);
            p(j).v  = vin(j);
        end
        p = p.normalize;
        p = reshape(p,size(pt));
        %p = normalize_plane(p);
    end %pointdirections()
    
    % ---------------------------------------------------
    function p = pointnormal(pt,nml)
        % Create plane from point (point3d) and normal (vector3d)
        assert( numel(pt) == numel(nml),'Number of elements must be the same')
        
        % Local anonymous functions:
        normalize_d = @(x)x./repmat( sqrt( sum( x.^2, 1 ) ), [size(x,1),1] );    
        
        % Plane normal, as 3-b-N double
        nmld = normalize_d([nml.uvw]);   % Convert to double & normalize
        [ndim,np] = size(nmld);                
        % To determine the in-plane u & v vectors, take uvw and add to the
        % smallest component to create a vector that is non-parallel with
        % nml:        
        [~,ls_rid] = min(abs(nmld),[],1);      % Get least significant component
        idx = sub2ind([ndim,np],ls_rid,(1:np));% Convert rowid to index
        vdash = nmld;                   % Initialise temporary vector
        vdash(idx) = vdash(idx) + 1;    % Select least significant comp. & add 1
        
        % Now the cross of this vector with nml will produce a
        % mutually-normal vector, which will be our first in-plane vector:
        uhat = normalize_d( cross(nmld,vdash) );
        % Then cross with nml to get the second in-plane vector:
        vhat = cross(uhat,nmld);
        % Convert to vector3d objects:
        siz = size(pt);
        uhat = reshape( vector3d(uhat), siz );
        vhat = reshape( vector3d(vhat), siz );
        % Plot check:
        %{
        figure, hold on
        quiver3(0,0,0,nmld(1,1),nmld(2,1),nmld(3,1),'k')
        quiver3(0,0,0,udash(1,1),udash(2,1),udash(3,1),'k--')
        quiver3(0,0,0,uhat(1,1),uhat(2,1),uhat(3,1),'b')
        quiver3(0,0,0,vhat(1,1),vhat(2,1),vhat(3,1),'r')
        axis equal, axis tight, grid on, view(3)
        %}
        % Now create plane:
        p = plane3d.pointdirections(pt,vhat,uhat);
    end %pointnormal()
    
    % ---------------------------------------------------
    function p = throughpoints(p1,p2,p3)        
        % Create direction vectors:
        v1 = vector3d( [p2.xyz] - [p1.xyz] );
        v2 = vector3d( [p3.xyz] - [p1.xyz] );
        p = plane3d.pointdirections(p1,v1,v2);
        p = reshape(p,size(p1));
    end %throughpoints()
    
end %methods (Static)
              
end %classdef


% ========================================================================
%   Other Helper functions


% ------------------------------------------------------------------------
function out = chop( in, tol )
% function out = chop( in, tol )
% Replace values that differ from an integer by <= tol by the integer
% Inputs:
%  in       input array
%  tol      tolerance, default = eps
% Output:
%  out      input array with integer replacements, if any
if ~exist( 'tol', 'var' ) || isempty( tol )
    tol = eps;
end
out = in;
rin = round( in );
lx  = abs( rin - in ) <= tol;
out(lx) = rin(lx);
end %chop()


% ------------------------------------------------------------------------
function varargout = drawplane(p,varargin)
% Draw plane on current or specified axis
varargout = {};
% Find the axes:
if isnumeric(varargin{1}) && ishandle(varargin{1}) && strcmpi(get(varargin{1},'Type'),'axes')
    ax = varargin{1};
    varargin = varargin(2:end);
else
    idx = find( strcmpi(varargin,'parent') );
    if ~isempty(idx)
        [~,ax] = varargin{[idx,idx+1]};
        varargin([idx,idx+1]) = [];
    else
        ax = gca;
    end
end

% Plot limits
lim = get(ax, 'xlim');
xmin = lim(1);
xmax = lim(2);
lim = get(ax, 'ylim');
ymin = lim(1);
ymax = lim(2);
lim = get(ax, 'zlim');
zmin = lim(1);
zmax = lim(2);

% Lines corresponding to cube edges
lineX00 = line3d('point-vector', [xmin ymin zmin], [1 0 0]);
lineX01 = line3d('point-vector', [xmin ymin zmax], [1 0 0]);
lineX10 = line3d('point-vector', [xmin ymax zmin], [1 0 0]);
lineX11 = line3d('point-vector', [xmin ymax zmax], [1 0 0]);

lineY00 = line3d('point-vector', [xmin ymin zmin], [0 1 0]);
lineY01 = line3d('point-vector', [xmin ymin zmax], [0 1 0]);
lineY10 = line3d('point-vector', [xmax ymin zmin], [0 1 0]);
lineY11 = line3d('point-vector', [xmax ymin zmax], [0 1 0]);

lineZ00 = line3d('point-vector', [xmin ymin zmin], [0 0 1]);
lineZ01 = line3d('point-vector', [xmin ymax zmin], [0 0 1]);
lineZ10 = line3d('point-vector', [xmax ymin zmin], [0 0 1]);
lineZ11 = line3d('point-vector', [xmax ymax zmin], [0 0 1]);

% Compute intersection points with each plane
piX00 = p.intersect(lineX00);
piX01 = p.intersect(lineX01);
piX10 = p.intersect(lineX10);
piX11 = p.intersect(lineX11);
piY00 = p.intersect(lineY00);
piY01 = p.intersect(lineY01);
piY10 = p.intersect(lineY10);
piY11 = p.intersect(lineY11);
piZ00 = p.intersect(lineZ00);
piZ01 = p.intersect(lineZ01);
piZ10 = p.intersect(lineZ10);
piZ11 = p.intersect(lineZ11);

% Concatenate points into one array
points = [...
    piX00;piX01;piX10;piX11; ...
    piY00;piY01;piY10;piY11; ...
    piZ00;piZ01;piZ10;piZ11];

% check validity: keep only points inside window
ac = 1e-14;
vx = points.x>=xmin-ac & points.x<=xmax+ac;
vy = points.y>=ymin-ac & points.y<=ymax+ac;
vz = points.z>=zmin-ac & points.z<=zmax+ac;
valid = vx & vy & vz;

% If there is no intersection point, exit.
if all(valid==0)
    disp('plane is outside the drawing window');
    h = [];
    
else
    % Intersection points
    pts = unique(points(valid, :));
    
    % The two coordinate lines of the plane
    d1 = line3d(p.p0,p.u);
    d2 = line3d(p.p0,p.v);
    
    % Position of intersection points in plane coordinates
    u1 = d1.project(pts);
    u2 = d2.project(pts);
    
    % Reorder vertices in the correct order
    ind = convhull(u1, u2);
    ind = ind(1:end-1);
    
    % Draw the patch
    if ~isempty(varargin)
        opts = varargin;
    else
        opts = {'m'};
    end
    h = patch(pts(ind).x, pts(ind).y, pts(ind).z, 'Parent',ax, opts{:});
end

% Return handle to plane if needed
if nargout>0
    varargout{1}=h;
end

end %drawplane()
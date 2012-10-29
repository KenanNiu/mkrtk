classdef vector3d
%VECTOR class.    
%
% Constructors:
%   v = vector3d()
%
% Methods:
%   angle(v1,v2)    % Angle between two vectors
%   cross(v1,v2)    % Cross product
%   dot(v1,v2)      % Dot product
%   double(v)       % Convert to double
%   eq(v1,v2)       % Test vector equality (==)
%   equiv(v1,v2)    % Test if vectors
%   neq(v1,v2)      % Not equal to (~=)
%   norm(v)         % Norm of vector
%   normalize(v)    % Normalize vector
%   plus(v1,v2)     % Vector addition
%   times(v,s)      % Vector scaling, multiplication by a scalar
%   u
%   v
%   w
%
%   See also POINT3D, LINE3D


properties (SetAccess = protected)
    uvw = [0 0 0]';
end % properties

% properties (Dependent = true, Hidden = true, SetAccess = private)
%     u
%     v
%     w
% end

methods % ordinary methods
    function v = vector3d( uvw ) % (constructor)
        switch nargin
            case 0  % nargin == 0
                %default uvw
            case 1  % nargin == 1
                siz = size(uvw);
                [arg3, dim3, perm] = finddim( uvw, 3 );
                if dim3 == 1  || dim3 == 2
                    siz(dim3) = 1;
                    nel       = prod( siz );
                    for iel = nel : -1 : 1
                        v(iel).uvw = arg3(:,iel);
                    end
                    v = reshape(  v, siz  );
                    v = ipermute( v, perm );
                    v = squeeze(v);   % Drop leading singleton dimension
                else
                    es = ['Incorrectly sized input. The three (x,y,z) ',...
                        'components should be on the first or second dimension.'];
                    error(es)
                end
                
            otherwise
                error('Incorrect number of input arguments')
        end %switch nargin
    end %vector3d()
    
    
    %=====================================================================
    % Ordinary Methods
    
    
    %----------------------------------------------------
    function theta = angle(v1,obj)
        switch class(obj)
            case 'vector3d'
                % Angle between two vectors.
                % http://www.mathworks.com/matlabcentral/newsreader/view_thread/151925#381952
                v2 = obj;
                a = [v1.uvw];
                b = [v2.uvw];
                theta = atan2( norm( cross(a,b) ), dot(a,b) );
                if numel(v1) > numel(v2)
                    siz = size(v1);
                else
                    siz = size(v2);
                end
                theta = reshape( theta, siz );
            case 'plane'
                % Minimal angle between vector and plane
                error('not yet implemented')
        end %switch
    end %angle()
    
    %----------------------------------------------------
    function v = cross(v1,v2)
        %CROSS Vector cross product
        %       V = V1 x V2
        %  V1 and V1 are N-by-M-by-... 
        %  If either V1 or V2 is 1-by-1, the result will be
        %  the size of the other input and will be the 
        %  cross product of this vector with each element.
        
        % Account for different input sizes:
        a = [v1.uvw];
        b = [v2.uvw];
        if numel(v1) > numel(v2)
            sizv = size(v1);
        else
            sizv = size(v2);
        end
        % Standard method: cross( [v1.uvw], [v2.uvw] );
        % Instead use faster method using bsxfun 
        c = cross(a,b); %helper function; not this class method
        v = reshape( vector3d( c ), sizv );
    end %cross()
    
    %----------------------------------------------------
    function d = dot(v1,v2)
        %DOT Scalar (dot) product
        %       D = V1 . V2
        %  V1 and V1 are N-by-M-by-... 
        %  If either V1 or V2 is 1-by-1, the result will be
        %  the size of the other input and will be the 
        %  dot product of this vector with each element.
        
        % Account for different input sizes:
        a = [v1.uvw];
        b = [v2.uvw];
        if numel(v1) > numel(v2)
            sizv = size(v1);
        else
            sizv = size(v2);
        end
        d = dot(a,b);
        d = reshape(d,sizv);
    end %dot()
    
    %----------------------------------------------------
    function d = double(v)
        %DOUBLE Convert to matrix format
        % If V is N-by-M-by..., then D will be 3-by-N-by-M-by...
        d = [v.uvw];
        d = reshape(d,[3,size(v)]);
    end %double()
    
    %----------------------------------------------------
    function tf = eq(v1,v2)
        if strcmpi( class(v1), class(v2) )
            si1 = size( v1 );
            si2 = size( v2 );
            ne1 = prod( si1 );
            ne2 = prod( si2 );
            v1_xyz = chop([v1.uvw]);
            v2_xyz = chop([v2.uvw]);
            if (ne1 == 0) || (ne2 == 0)
                tf   = logical([]);
                return;
            elseif ne1 == 1
                tf   = repmat( v1_xyz, 1, ne2 ) == v2_xyz;
                siz = si2;
            elseif ne2 == 1
                tf   = v1_xyz == repmat( v2_xyz, 1, ne1 );
                siz = si1;
            elseif isequal( si1, si2 )
                tf   = v1_xyz == v2_xyz;
                siz = si1;
            else
                error( 'Matrix dimensions must agree' );
            end
            tf   = reshape( all( tf, 1 ), siz );
        else
            tf = false;
        end
    end %eq()
    
    %----------------------------------------------------
    function tf = equiv(v1,v2)
        % EQUIV Test vector equivalence - are they parallel?
        if strcmpi( class(v1), class(v2) )
            tf = eq( v1.normalize, v2.normalize );
        else
            tf = false;
        end
    end %equiv()
    
%     %----------------------------------------------------
%     function tf = isparallel(v1,v2)
%         % Test if vectors are parallel - alias for equiv()
%         tf = equiv(v1,v2);
%     end %isparallel()
    
    %----------------------------------------------------
    function v = mtimes(v,a)
        if isnumeric(a) && numel(a)==1
            v = times(v,a);
        elseif ~isnumeric(a)
            error('Matrix multiplication not defined for class ''%s''',class(a))
        else
            error('Matrix multiplication not yet defined.')
        end
    end %mtimes()
    
    %----------------------------------------------------
    function tf = ne(v1,v2)
        tf = ~eq( v1, v2 );
    end %ne()
    
    %----------------------------------------------------
    function n = norm(v)
        %NORM Vector magnitude (norm)        
        n = norm([v.uvw]); % Calls helper function by the same name
        n = reshape(n,size(v));
    end %norm()
    
    %----------------------------------------------------
    function vhat = normalize(v)
        % Enforce magnitued == 1
        n = v.norm;
        if all(n == 1)
            vhat = v;   % short circuit if already done
        else
            N = repmat(n(:)',[3,1]);
            vhat = vector3d( [v.uvw]./N );
            vhat = reshape( vhat, size(v) );
        end
    end %normalize()

    %----------------------------------------------------
    function v = plus(v1,v2)
        if strcmpi( class(v1), class(v2) )
            si1 = size( v1 );
            si2 = size( v2 );
            ne1 = prod( si1 );
            ne2 = prod( si2 );
            %if (ne1 == 0) || (ne2 == 0)
            %    v = v1;
            if ne1 == 1
                v = [v2.uvw] + repmat( [v1.uvw], 1, ne2 );
                siz = si2;
            elseif ne2 == 1
                v = [v1.uvw] + repmat( [v2.uvw], 1, ne1 );
                siz = si1;
            elseif isequal( si1, si2 )
                v = [v1.uvw] + [v2.uvw];
                siz = si1;
            else
                error( 'Matrix dimensions must agree' );
            end
            v = reshape( vector3d(v), siz );
        else
            error(['Only vector3d objects can be added to vector3d objects.  ',...
                'Adding a ''%s'' to a vector3d object is not allowed'],class(v2))
        end
    end %plus()
        
    %----------------------------------------------------
    function vs = times(v,s)
        if numel(s) == 1 && isnumeric(s)
            vs = vector3d( [v.uvw].*s );
            vs = reshape( vs, size(v) );
        else
            error('Element-by-element multiplication not yet supported')
        end
    end %times()
    
    
    % ===================================================
    % The following methods act similarly to getter methods.
    % However, we do not use getter methods because they act 
    % on each individual object, rather than on an array of
    % objects.  We want to act on the array of objects and
    % return arrays as outputs.
    
    %----------------------------------------------------
    function u1 = u(vec)
        d = [vec.uvw];
        u1 = d(1,:);
        u1 = reshape(u1,size(vec));
    end %u()
    
    %----------------------------------------------------
    function v1 = v(vec)
        d = [vec.uvw];
        v1 = d(2,:);
        v1 = reshape(v1,size(vec));
    end %v()
    
    %----------------------------------------------------
    function w1 = w(vec)
        d = [vec.uvw];
        w1 = d(3,:);
        w1 = reshape(w1,size(vec));
    end %w()
    
end %methods



methods (Static)
end %methods (Static)

end %classdef vector3d

% ========================================================================
%   Helper functions


% ------------------------------------------------------------------------
function [aout, dim, perm] = finddim( ain, len )
% Find first dimension in ain of length len, permute ain to make it first
% Inputs:
%  ain(s1,s2,...)   data array, size = [s1, s2, ...]
%  len              length sought, e.g. s2 == len
%                   if len < 0, then find first dimension >= |len|
% Outputs:
%  aout(s2,...,s1)  data array, permuted so first dimension is length len
%  dim              dimension number of length len, 0 if ain has none
%  perm             permutation order (for permute and ipermute) of aout,
%                   e.g. [2, ..., 1]
% Notes: if no dimension has length len, aout = ain, dim = 0, perm = 1:ndm
%        ain = ipermute( aout, perm )
siz  = size( ain );
ndm  = length( siz );
if len < 0
    dim  = find( siz >= -len, 1, 'first' );
else
    dim  = find( siz == len, 1, 'first' );
end
if isempty( dim )
    dim  = 0;
end
if dim < 2
    aout = ain;
    perm = 1 : ndm;
else
    % Permute so that dim becomes the first dimension
    %perm = [ dim : ndm, 1 : dim-1 ]; % <- equivalent to shiftdim
    perm = [dim, 1:ndm]; 
    perm(dim+1) = [];
    aout = permute( ain, perm );
end
end %finddim()


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
function c = cross(a,b)
%CROSS_PRODUCT Vector cross product using bsxfun (faster that CROSS)
%  Calculate the vector cross product as follows:
%       C = A x B
%
%  A and B must be M-by-3.  Alternatively, one can be 1-by-3, in which case
%  the output will be the cross product of this with all the rows of the
%  other input.
%
%  WARNING!!  This function does not act on vector3d class objects, but on
%  float arrays. We have a class method by the same name, CROSS, for
%  objects. 
%  To call that method, use the a.cross(b) syntax on the object, instead of
%  the cross(a,b) syntax.
%
assert( ~isa(a,'vector3d') && ~isa(b,'vector3d'),...
    'Use a.cross(b) to calculate the cross product of two vector3d objects')

assert(size(a,1)==3,'A must be 3-by-N');
assert(size(b,1)==3,'B must be 3-by-N');
siza = size(a);
sizb = size(b);
if siza(2) == 1
    sizc = size(b);
elseif sizb(2) == 1
    sizc = size(a);
elseif siza(2) == sizb(2)
    sizc = size(a);
else
    error('Matrix dimensions must agree')
end
c = zeros(sizc);
c(:) = bsxfun(@times, a([2 3 1],:,:), b([3 1 2],:,:)) - ...
       bsxfun(@times, b([2 3 1],:,:), a([3 1 2],:,:));
end %cross_product()


% ------------------------------------------------------------------------
function d = dot(a,b)
%DOT Vector dot (scalar) product
%
%  WARNING!!  This function does not act on vector3d class objects, but on
%  float arrays.  We have a class method by the same name, DOT, for
%  objects. 
%  To call that method, use the a.dot(b) syntax on the object, instead of the
%  dot(a,b) syntax.
%
assert( ~isa(a,'vector3d') && ~isa(b,'vector3d'),...
    'Use a.dot(b) to calculate the dot product of two vector3d objects')
assert(size(a,1)==3,'A must be 3-by-N');
assert(size(b,1)==3,'B must be 3-by-N');

d = sum(bsxfun(@times, a, b),1);

end %dot()


% ------------------------------------------------------------------------
function n = norm(x)
%NORM Vector norm
%
%  WARNING!!  This function does not act on vector3d class objects, but on
%  float arrays.  We have a class method by the same name, NORM, for
%  objects. 
%  To call that method, use the v.norm syntax on the object, instead of the
%  norm(x) syntax.
%
assert(~isa(x,'vector3d'),'Use v.norm to calculate the norm of a vector3d object')

assert(size(x,1) == 3,'X must be 3-by-N')
n = sqrt(sum(x.*x,1));
end %norm()

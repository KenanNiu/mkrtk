classdef point3d
%POINT3D class.
%
% Constructors:
%   v = point3d()
%
% Methods:
%   angle(p1,p2,p3)     % Angle between three points
%   double(p)           % Convert to double
%   eq(p1,p2)           % Test vector equality (==)
%   neq(p1,p2)          % Not equal to (~=)
%   plus(p1,xyz)          % Vector addition; essentially amounts to a translation
%   rotate(p,???)       % Inputs?
%   times(p,s)          % Multiplication by a scalar
%   translate(p,d)      % Translate
%   x
%   y
%   z
%
% See also VECTOR3D, LINE3D


properties (SetAccess = protected)
    xyz = [0 0 0]';
end % properties

% properties (Dependent = true, Hidden = true, SetAccess = private)
%     x
%     y
%     z
% end

methods % ordinary methods
    function p = point3d( xyz ) % (constructor)
        switch nargin
            case 0  % nargin == 0
                %default xyz
            case 1  % nargin == 1
                siz = size(xyz);
                [arg3, dim3, perm] = finddim( xyz, 3 );
                if dim3 == 1  || dim3 == 2
                    siz(dim3) = 1;
                    nel       = prod( siz );
                    for iel = nel : -1 : 1
                        p(iel).xyz = arg3(:,iel);
                    end
                    p = reshape(  p, siz  );
                    p = ipermute( p, perm );
                    p = squeeze(p);   % Drop leading singleton dimension
                else
                    es = ['Incorrectly sized input. The three (x,y,z) ',...
                        'components should be on the first or second dimension.'];
                    error(es)
                end
                
            otherwise
                error('Incorrect number of input arguments')
        end %switch nargin
    end %point3d()
    
    
    %=====================================================================
    % Ordinary Methods
    
    
    %----------------------------------------------------
    function theta = angle(p1,p2,p3)
        % ANGLE Calculate angle between 3 points
        % The angle is p1-p2-p3 where p2 is the apex
        if strcmpi( class(p1), class(p2) ) && strcmpi( class(p1), class(p3) )
            siz1 = size(p1);
            siz2 = size(p2);
            siz3 = size(p3);
            p1_xyz = [p1.xyz];
            p2_xyz = [p2.xyz];
            p3_xyz = [p3.xyz];
            if isequal( siz1, siz2 ) && isequal( siz1, siz3 )
                % All same size
                siz = siz1;
            elseif numel(p1)==1 && isequal( siz2, siz3 )
                % #1 for repeat:
                siz = siz2;
                p1_xyz = repmat(p1.xyz,[1,prod(siz)]);
            elseif numel(p2)==1 && isequal( siz1, siz3 )
                % #2 for repeat 
                siz = siz1;
                p2_xyz = repmat(p2.xyz,[1,prod(siz)]);
            elseif numel(p3)==1 && isequal( siz1, siz2 ) 
                % #3 for repeat
                siz = siz1;
                p3_xyz = repmat(p3.xyz,[1,prod(siz)]);
            else
                error('Unhandled input dimensions')
            end
            v1 = p1_xyz - p2_xyz;                   %\_ Create direction vectors
            v2 = p3_xyz - p2_xyz;                   %/
            norm = @(v)sqrt( sum( v.^2,1 ) );       %\
            v1 = v1./repmat( norm(v1), [3,1] );     % |- Normalize vectors
            v2 = v2./repmat( norm(v2), [3,1] );     %/
            theta = acos( dot( v1, v2, 1 ) );       %- Calculate angle
            theta = reshape(theta, siz);            %- Restore shape
        else
            error('Inputs must be point3d objects')
        end
    end %angle()
    

    
    %----------------------------------------------------
    function d = double(p)
        %DOUBLE Convert to matrix format
        % If V is N-by-M-by..., then D will be 3-by-N-by-M-by...
        d = [p.xyz];
        d = reshape(d,[3,size(p)]);
    end %double()
    
    %----------------------------------------------------
    function tf = eq(p1,p2)
        if strcmpi( class(p1), class(p2) )
            si1 = size( p1 );
            si2 = size( p2 );
            ne1 = prod( si1 );
            ne2 = prod( si2 );
            p1_xyz = chop([p1.xyz]);
            p2_xyz = chop([p2.xyz]);
            if (ne1 == 0) || (ne2 == 0)
                tf   = logical([]);
                return;
            elseif ne1 == 1
                tf   = repmat( p1_xyz, 1, ne2 ) == p2_xyz;
                siz = si2;
            elseif ne2 == 1
                tf   = p1_xyz == repmat( p2_xyz, 1, ne1 );
                siz = si1;
            elseif isequal( si1, si2 )
                tf   = p1_xyz == p2_xyz;
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
    function p = mtimes(p,a)
        if isnumeric(a) && numel(a)==1
            p = times(p,a);
        elseif ~isnumeric(a)
            error('Matrix multiplication not defined for class ''%s''',class(a))
        else
            error('Matrix multiplication not yet defined.')
        end
    end %mtimes()
    
    %----------------------------------------------------
    function tf = ne(p1,p2)
        tf = ~eq( p1, p2 );
    end %ne()

    %----------------------------------------------------
    function ps = plus(p,s)
        if numel(s) == 3 && isnumeric(s)
            ps = point3d( [p.xyz] + repmat( s(:), [1, numel(p)] ) );
            ps = reshape( ps, size(p) );
        else
            error('Only [x,y,z] displacements can be added to point3d objects')
        end
    end %plus()
        
    %----------------------------------------------------
    function ps = times(p,s)
        if numel(s) == 1 && isnumeric(s)
            ps = point3d( [p.xyz].*s );
        elseif numel(s) == 3 && isnumeric(s)
            ps = point3d( [p.xyz].* repmat( s(:), [1, numel(p)] ) );
        else
            error('Element-by-element multiplication not yet supported')
        end
        ps = reshape( ps, size(p) );
    end %times()
    
    %----------------------------------------------------
    function ps = translate(p,s)
        % Translate: alias for plus()
        ps = plus(p,s);
    end %translate()
    
    %----------------------------------------------------
    function pu = unique(p)
        pu = point3d(unique([p.xyz]','rows')');
    end %unique()
    
    % ===================================================
    % The following methods act similarly to getter methods.
    % However, we do not use getter methods because they act 
    % on each individual object, rather than on an array of
    % objects.  We want to act on the array of objects and
    % return arrays as outputs.
    
    %----------------------------------------------------
    function x1 = x(p)
        d = [p.xyz];
        x1 = d(1,:);
        x1 = reshape(x1,size(p));
    end %x()
    
    %----------------------------------------------------
    function y1 = y(p)
        d = [p.xyz];
        y1 = d(2,:);
        y1 = reshape(y1,size(p));
    end %y()
    
    %----------------------------------------------------
    function z1 = z(p)
        d = [p.xyz];
        z1 = d(3,:);
        z1 = reshape(z1,size(p));
    end %z()
    
end %methods



methods (Static)
end %methods (Static)

end %classdef vector

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

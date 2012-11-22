classdef quaternion
%QUATERNION class.  
%
%
% Properties (SetAccess = protected):
%  e(4,1)   components, basis [1; i; j; k]: e(1) + i*e(2) + j*e(3) + k*e(4)
%           i*j=k, j*i=-k, j*k=i, k*j=-i, k*i=j, i*k=-j, i*i = j*j = k*k = -1
%
%
% Constructors:
%   q = quaternion              
%   q = quaternion(X)           
%   q = quaternion(w,x,y,z)     
%
%
% This file is derived from the quaternion class by Mark Tincknell, which
% can be found here:
% http://www.mathworks.com/matlabcentral/fileexchange/33341-quaternion-m 
%
% Last cross-checked version: 13 Nov 2012
%
% Changes include, but are not limited to
%   - reduced method set, using only those required for this application
%   - some data in different formats or order
%   - some extra methods where required
%   - some corrections or changes in convention
%
% Revised/changed by Joshua Martin, Nov-2012


properties (SetAccess = protected)
    e   = zeros(4,1);
end % properties


methods
    %---------------------------------------------------
    % Array constructors
    function q = quaternion( varargin ) % (constructor)
        switch nargin

            case 0  % nargin == 0
                q.e = zeros(4,1);
                return;

            case 1  % nargin == 1
                siz = size( varargin{1} );
                nel = prod( siz );
                if nel == 0
                    q	= quaternion.empty;
                    return;
                elseif nel == 1
                    q.e = chop( [varargin{1}; zeros(3,1)] );
                    return
                elseif isa( varargin{1}, 'quaternion' )
                    q   = varargin{1};
                    return;
                elseif (nel == 1) || ~isreal( varargin{1}(:) )
                    for iel = nel : -1 : 1
                        q(iel).e = chop( [real(varargin{1}(iel)); ...
                                          imag(varargin{1}(iel)); ...
                                          0; ...
                                          0] );
                    end
                    q   = reshape( q, siz );
                    return;
                end
                [arg4, dim4, perm4] = finddim( varargin{1}, 4 );
                if dim4 > 0
                    siz(dim4)   = 1;
                    nel         = prod( siz );
                    if dim4 > 1
                        perm    = perm4;
                    end
                    for iel = nel : -1 : 1
                        q(iel).e = chop( arg4(:,iel) );
                    end
                else
                    % Imaginary components
                    [arg3, dim3, perm3] = finddim( varargin{1}, 3 );
                    if dim3 > 0
                        siz(dim3)   = 1;
                        nel         = prod( siz );
                        if dim3 > 1
                            perm    = perm3;
                        end
                        for iel = nel : -1 : 1
                            q(iel).e = chop( [0; arg3(:,iel)] );
                        end
                    else
                        error( 'Invalid input' );
                    end
                end
                
            case 2  % nargin == 2
                % real-imaginary only (no j or k) inputs
                varargin    = conformsize( varargin{:} );
                siz         = size( varargin{1} );
                nel         = prod( siz );
                for iel = nel : -1 : 1
                    q(iel).e = chop( [varargin{1}(iel); ...
                                      varargin{2}(iel); ...
                                      0;
                                      0] );
                end
                
            case 3  % nargin == 3
                % vector inputs (no real, only i, j, k)
                varargin    = conformsize( varargin{:} );
                siz         = size( varargin{1} );
                nel         = prod( siz );
                for iel = nel : -1 : 1
                    q(iel).e = chop( [0; ...
                                      varargin{1}(iel); ...
                                      varargin{2}(iel); ...
                                      varargin{3}(iel)] );
                end

            otherwise   % nargin >= 4
                varargin    = conformsize( varargin{:} );
                siz         = size( varargin{1} );
                nel         = prod( siz );
                for iel = nel : -1 : 1
                    q(iel).e = chop( [varargin{1}(iel); ...
                                      varargin{2}(iel); ...
                                      varargin{3}(iel); ...
                                      varargin{4}(iel)] );
                end
        end % switch nargin

        q   = reshape( q, siz );
        if exist( 'perm', 'var' )
            q   = ipermute( q, perm );
        end
    end % quaternion (constructor)
    
    
    %=====================================================================
    % Ordinary Methods
    
    %---------------------------------------------------
    function n = abs( q )
        % Length of quaternion
        n = q.modulus;
    end %abs

	%---------------------------------------------------
    function q3 = bsxfun( func, q1, q2 )
        % Binary Singleton Expansion for quaternion arrays. Apply the element by
        % element binary operation specified by the function handle func to arrays
        % q1 and q2. All dimensions of q1 and q2 must either agree or be length 1.
        % Inputs:
        %  func     function handle (e.g. @plus) of quaternion function or operator
        %  q1(n1)   quaternion array
        %  q2(n2)   quaternion array
        % Output:
        %  q3(n3)   quaternion array of function or operator outputs
        %           size(q3) = max( size(q1), size(q2) )
        [q1,q2] = qcheck('assign',q1,q2);
        s1  = size( q1 );
        s2  = size( q2 );
        nd1 = length( s1 );
        nd2 = length( s2 );
        s1  = [s1, ones(1,nd2-nd1)];
        s2  = [s2, ones(1,nd1-nd2)];
        if ~all( (s1 == s2) | (s1 == 1) | (s2 == 1) )
            error( 'Non-singleton dimensions of q1 and q2 must match each other' );
        end
        c1  = num2cell( s1 );
        c2  = num2cell( s2 );
        s3  = max( s1, s2 );
        nd3 = length( s3 );
        n3  = prod( s3 );
        q3  = quaternion.nan( s3 );
        for i3 = 1 : n3
            [ix3{1:nd3}] = ind2sub( s3, i3 ); 
            ix1     = cellfun( @min, ix3, c1, 'UniformOutput', false );
            ix2     = cellfun( @min, ix3, c2, 'UniformOutput', false );
            q3(i3)  = func( q1(ix1{:}), q2(ix2{:}) );
        end
    end % bsxfun

	%---------------------------------------------------
    function qc = conj( q )
        % Conjugate: qc = [q(1) -q(2) -q(3) -q(4)]'
        d   = double( q );
        d(2:4,:) = -d(2:4,:);
        qc  = reshape( quaternion( d ), size( q ));
    end %conj
    
    %---------------------------------------------------
    function qt = ctranspose( q )
        % Conjugate transpose
        qt  = transpose( q.conj );
    end % ctranspose
    
    %---------------------------------------------------
    function qd = diff( q, ord, dim )
        % Difference between quaternions (array difference)
        % quaternion array difference, ord is the order of difference (default = 1)
        % dim defaults to first dimension of length > 1
        if isempty( q )
            qd  = q;
            return;
        end
        if ~exist( 'ord', 'var' ) || isempty( ord )
            ord = 1;
        end
        if ord <= 0
            qd  = q;
            return;
        end
        if ~exist( 'dim', 'var' )
            [q, dim, perm]  = finddim( q, -2 );
        elseif dim > 1
            ndm  = ndims( q );
            perm = [ dim : ndm, 1 : dim-1 ];
            q    = permute( q, perm );
        end
        siz = size( q );
        if siz(1) <= 1
            qd  = quaternion.empty;
            return;
        end
        for is = siz(1)-1 : -1 : 1
            qd(is,:) = q(is+1,:) - q(is,:);
        end
        ord = ord - 1;
        if ord > 0
            qd  = diff( qd, ord, 1 );
        end
        if dim > 1
            qd  = ipermute( qd, perm );
        end
    end % diff
    
    %---------------------------------------------------
    function display( q )
        siz = size( q );
        nel = [1 cumprod( siz )];
        ndm = length( siz );
        disp(' ')
        fprintf('%s = \n\n',inputname(1));
        sstr = sprintf('%dx',siz);  sstr(end) = [];
        fprintf('  %s <a href="matlab:help quaternion">quaternion</a>',sstr);
        fprintf(' with <a href="matlab:methods(''quaternion'')">methods</a>\n\n');
        for iel = 1 : nel(end)
            if nel(end) == 1
                sub = '';
            else
                sub = ')\t';
                jel = iel - 1;
                for idm = ndm : -1 : 1
                    idx = floor( jel / nel(idm) ) + 1;
                    sub = [',' int2str(idx) sub]; %#ok<AGROW>
                    jel = rem( jel, nel(idm) );
                end
                sub(1)  = '(';
            end
            fprintf( ['%s' sub ' = %7.4g + %7.4g i + %7.4g j + %7.4g k\n'], ...
                inputname(1), q(iel).e )
        end
        fprintf('\n');
        %fprintf('\n  <a href="matlab:methods(''quaternion'')">Methods</a>\n\n');
        
    end
    
    %---------------------------------------------------
    function d = double( q )
        siz = size( q );
        d   = reshape( [q.e], [4 siz] );
        d   = chop( d );
    end % double()
    
    %---------------------------------------------------
    function l = eq( q1, q2 )
        if isa( q1, 'quaternion' ) && isa( q2, 'quaternion' )
            si1 = size( q1 );
            si2 = size( q2 );
            ne1 = prod( si1 );
            ne2 = prod( si2 );
            if (ne1 == 0) || (ne2 == 0)
                l   = logical([]);
                return;
            elseif ne1 == 1
                siz = si2;
            elseif ne2 == 1
                siz = si1;
            elseif isequal( si1, si2 )
                siz = si1;
            else
                error( 'Matrix dimensions must agree' );
            end
        else
            error( 'Inputs must be quaternions' );
        end
        l   = bsxfun( @eq, [q1.e], [q2.e] );
        l   = reshape( all( l, 1 ), siz );
    end % eq()
    
    %---------------------------------------------------
    function l = equiv( q1, q2, tol )
        % Quaternion rotational equivalence, within tolerance tol, 
        %   l = (q1 == q2) | (q1 == -q2)
        % Optional argument tol (default = eps) sets tolererance for
        % difference from exact equality
        if ~exist( 'tol', 'var' )
            tol = eps;
        end
        si1 = size( q1 );
        si2 = size( q2 );
        ne1 = prod( si1 );
        ne2 = prod( si2 );
        if (ne1 == 0) || (ne2 == 0)
            l   = logical([]);
            return;
        elseif ne1 == 1
            siz = si2;
        elseif ne2 == 1
            siz = si1;
        elseif isequal( si1, si2 )
            siz = si1;
        else
            error( 'Matrix dimensions must agree' );
        end
        dm  = chop( bsxfun( @minus, [q1.e], [q2.e] ), tol );
        dp  = chop( bsxfun( @plus,  [q1.e], [q2.e] ), tol );
        l   = all( (dm == 0) | (dp == 0), 1 );
        l   = reshape( l, siz );
    end % equiv()
    
    %---------------------------------------------------
    function qi = inverse( q )
        % Quaternion inverse:
        %   qi = conj(q)/norm(q)^2, q*qi = qi*q = 1 for q ~= 0
        if isempty( q )
            qi  = q;
            return;
        end
        d   = double( q );
        d(2:4,:) = -d(2:4,:);
        n2  = repmat( sum( d.^2, 1 ), 4, ones( 1, ndims( d ) - 1 ));
        ne0 = n2 ~= 0;
        di  = Inf( size( d ));
        di(ne0)  = d(ne0) ./ n2(ne0);
        qi  = reshape( quaternion( di ), size( q ));
    end % inverse()
    
    %---------------------------------------------------
    function l = isfinite( q )
        % Check if all quatenrion elements are finite
        l   = reshape( all( isfinite( [q.e] ), 1 ), size( q ));
    end % isfinite()

    %---------------------------------------------------
    function l = isinf( q )
        % Check if quaternion has any Inf components
        l   = reshape( any( isinf( [q.e] ), 1 ), size( q ));
    end % isinf()

    %---------------------------------------------------
    function l = isnan( q )
        % Check if quaternion has any NaN components
        l = reshape(any(isnan([q.e])),size(q));
    end % isnan()
    
    %---------------------------------------------------
    function q3 = ldivide( q1, q2 )
        if ~isa( q1, 'quaternion' )
            q1  = quaternion( q1, 0, 0, 0 );
        end
        if ~isa( q2, 'quaternion' )
            q2  = quaternion( q2, 0, 0, 0 );
        end
        si1 = size( q1 );
        si2 = size( q2 );
        ne1 = prod( si1 );
        ne2 = prod( si2 );
        if (ne1 == 0) || (ne2 == 0)
            q3  = quaternion.empty;
            return;
        elseif ~isequal( si1, si2 ) && (ne1 ~= 1) && (ne2 ~= 1)
            error( 'Matrix dimensions must agree' );
        end
        if ne2 > ne1
            q3  = repmat( quaternion, si2 );
        else
            q3  = repmat( quaternion, si1 );
        end
        for iel = 1 : max( ne1, ne2 )
            q3(iel) = product( q1(min(iel,ne1)).inverse, ...
                               q2(min(iel,ne2)) );
        end
    end % ldivide()
    
    %---------------------------------------------------
    function q3 = minus( q1, q2 )
        % Quaternion subtraction: q3 = q1 - q2
        if ~isa( q1, 'quaternion' )
            q1  = quaternion( q1, 0, 0, 0 );
        end
        if ~isa( q2, 'quaternion' )
            q2  = quaternion( q2, 0, 0, 0 );
        end
        si1 = size( q1 );
        si2 = size( q2 );
        ne1 = prod( si1 );
        ne2 = prod( si2 );
        if (ne1 == 0) || (ne2 == 0)
            q3  = quaternion.empty;
            return;
        elseif ne1 == 1
            siz = si2;
        elseif ne2 == 1
            siz = si1;
        elseif isequal( si1, si2 )
            siz = si1;
        else
            error( 'Matrix dimensions must agree' );
        end
        d3  = bsxfun( @minus, [q1.e], [q2.e] );
        q3  = quaternion( d3 );
        q3  = reshape( q3, siz );
    end % minus()
    
    %---------------------------------------------------
    function m = modulus( q )
        % Modulus of quaternion
        m = q.norm;
    end % modulus()
    
    %---------------------------------------------------
    function q3 = mtimes( q1, q2 )
        % Matrix quaternion product of 2-D conformable quaternion matrices
        % q1 and q2:
        %   q3 = mtimes(q1,q2), or
        %   q3 = q1*q2
        [q1,q2] = qcheck('assign',q1,q2);
        si1 = size( q1 );
        si2 = size( q2 );
        ne1 = prod( si1 );
        ne2 = prod( si2 );
        if (ne1 == 1) || (ne2 == 1)
            q3  = times( q1, q2 );
            return;
        end
        if (length( si1 ) ~= 2) || (length( si2 ) ~= 2)
            error( 'Input arguments must be 2-D' );
        end
        if si1(2) ~= si2(1)
            error( 'Inner matrix dimensions must agree' );
        end
        q3  = repmat( quaternion, [si1(1) si2(2)] );
        for i1 = 1 : si1(1)
            for i2 = 1 : si2(2)
                for i3 = 1 : si1(2)
                    q3(i1,i2) = q3(i1,i2) + product( q1(i1,i3), q2(i3,i2) );
                end
            end
        end
    end % mtimes()

    %---------------------------------------------------
    function l = ne( q1, q2 )
        % Not equal
        l   = ~eq( q1, q2 );
    end % ne()
        
    %---------------------------------------------------
    function n = norm( q )
        % Norm of quaternion
        % This can be defined as
        %       sqrt( w^2 + x^2 + y^2 + z^2 )
        %   or as
        %       w^2 + x^2 + y^2 + z^2 
        %   and can be using the conjugate:
        %       (q*qconj) == (qconj*q)
        %       
        % Here we use the first definition, ie, norm == modulus
        n = shiftdim( sqrt( sum( double( q ).^2, 1 )), 1 );
    end %norm()
    
    %---------------------------------------------------
    function [q, n] = normalize( q )
        % Normalize quaternion
        %   [q, n] = normalize( q )
        % q = quaternions with norm == 1 (unless q == 0), n = former norms
        siz = size( q );
        nel = prod( siz );
        ndm = length( siz );
        if nel == 0
            if nargout > 1
                n   = zeros( siz );
            end
            return;
        elseif nel > 1
            nel = [];
        end
        d   = double( q );
        n   = sqrt( sum( d.^2, 1 ));
        n4  = repmat( n, 4, nel );
        ne0 = (n4 ~= 0) & (n4 ~= 1);
        d(ne0)  = d(ne0) ./ n4(ne0);
        neg     = repmat( reshape( d(1,:) < 0, [1 siz] ), ...
                          [4, ones(1,ndm)] );
        d(neg)  = -d(neg);
        q       = reshape( quaternion( d ), siz );
        if nargout > 1
            n   = shiftdim( n, 1 );
        end
    end % normalize()
        
    %---------------------------------------------------
    function q3 = plus( q1, q2 )
        [q1,q2] = qcheck('assign',q1,q2);
        si1 = size( q1 );
        si2 = size( q2 );
        ne1 = prod( si1 );
        ne2 = prod( si2 );
        if (ne1 == 0) || (ne2 == 0)
            q3  = quaternion.empty;
            return;
        elseif ne1 == 1
            siz = si2;
        elseif ne2 == 1
            siz = si1;
        elseif isequal( si1, si2 )
            siz = si1;
        else
            error( 'Matrix dimensions must agree' );
        end
        d3  = bsxfun( @plus, [q1.e], [q2.e] );
        q3  = quaternion( d3 );
        q3  = reshape( q3, siz );
    end % plus()
    
    %---------------------------------------------------
    function qp = prod( q, dim )
        % Quaternion array product over dimension dim
        % dim defaults to first dimension of length > 1
        if isempty( q )
            qp  = q;
            return;
        end
        if ~exist( 'dim', 'var' )
            [q, dim, perm]  = finddim( q, -2 );
        elseif dim > 1
            ndm  = ndims( q );
            perm = [ dim : ndm, 1 : dim-1 ];
            q    = permute( q, perm );
        end
        siz = size( q );
        qp  = reshape( q(1,:), [1 siz(2:end)] );
        for is = 2 : siz(1)
            qp(1,:) = qp(1,:) .* q(is,:);
        end
        if dim > 1
            qp  = ipermute( qp, perm );
        end
    end % prod()
    
    %---------------------------------------------------
    function q3 = product( q1, q2 )
        % Quaternion product of scalar quaternions q1 and q2:
        %   q3 = product(q1,q2)
        [q1,q2] = qcheck('assign',q1,q2);
        if (numel( q1 ) ~= 1) || (numel( q2 ) ~= 1)
            error( 'product not defined for arrays, use mtimes or times' );
        end
        ee  = q1.e * q2.e.';
        eo  = [ee(1,1) - ee(2,2) - ee(3,3) - ee(4,4); ...
               ee(1,2) + ee(2,1) + ee(3,4) - ee(4,3); ...
               ee(1,3) - ee(2,4) + ee(3,1) + ee(4,2); ...
               ee(1,4) + ee(2,3) - ee(3,2) + ee(4,1)];
        q3  = quaternion( chop( eo ));
    end % product()
    
    %---------------------------------------------------
    function q3 = rdivide( q1, q2 )
        % This function does not produce the same result as the Aerospace
        % toolbox function QUATDIVIDE.  Should it?
        % Compare results of these two lines:
        %   >> (q(1)./q(2))
        %   >> quatdivide(squeeze([q(1).double])',squeeze([q(2).double])')
        [q1,q2] = qcheck('assign',q1,q2);
        si1 = size( q1 );
        si2 = size( q2 );
        ne1 = prod( si1 );
        ne2 = prod( si2 );
        if (ne1 == 0) || (ne2 == 0)
            q3  = quaternion.empty;
            return;
        elseif ~isequal( si1, si2 ) && (ne1 ~= 1) && (ne2 ~= 1)
            error( 'Matrix dimensions must agree' );
        end
        for iel = max( ne1, ne2 ) : -1 : 1
            q3(iel) = product( q1(min(iel,ne1)), ...
                               q2(min(iel,ne2)).inverse );
        end
        if ne2 > ne1
            q3  = reshape( q3, si2 );
        else
            q3  = reshape( q3, si1 );
        end
    end %rdivide()
    
    %---------------------------------------------------
    function qs = smooth( q , fun , dim)
        % This function smooths q along dimension dim with function handle
        % fun
        % Usage:
        %   qs = smooth(q, fun)
        %   qs = smooth(q, fun, dim)
        %
        % Example:
        %   n = 100;
        %   theta = linspac(0,2*pi,n);
        %   sfun  = @(y)smooth(theta,y,10/n,'rloess')
        %   qs = smooth(q, fun)         
        if numel(q) < 3
            error('Must have more than 3 quaternions')
        end
        if ~exist( 'dim', 'var' )
            [q, dim, perm]  = finddim( q, -2 );
        elseif dim > 1
            ndm  = ndims( q );
            perm = [ dim : ndm, 1 : dim-1 ];
            q    = permute( q, perm );
        end
        % Now reshape into 4-by-N-by-M double
        siz = size(q);
        q1 = reshape( q , siz(1),[] );
        d = q1.double;
        % Smooth components for each set
        ds = d;
        n = size(d,3);
        for j = 1:n
            for r = 1:4
                ds(r,:,j) = fun(d(r,:,j));
            end
        end
        % Un-shape:
        qs = reshape( quaternion(ds) , siz );
        if (size( qs, 1 ) == 1) 
            qs   = shiftdim( qs, 1 );
        end
        if dim > 1
            qs  = ipermute( qs, perm );
        end
    end %smooth()
    
    %---------------------------------------------------
    function q3 = times( q1, q2 )
        % Element-by-element quaternion multiplication:
        %   q3 = q1.*q2, or 
        %   q3 = times(q1,q2)
        %
        % This is the normal form for quaternion multiplication, 
        % and is equivalent to the Aerospace Toolbox command:
        %   quatmultiply(squeeze([q.double])',squeeze([q.double])')
        %
        [q1,q2] = qcheck('assign',q1,q2);
        si1 = size( q1 );
        si2 = size( q2 );
        ne1 = prod( si1 );
        ne2 = prod( si2 );
        if (ne1 == 0) || (ne2 == 0)
            q3  = quaternion.empty;
            return;
        elseif ~isequal( si1, si2 ) && (ne1 ~= 1) && (ne2 ~= 1)
            error( 'Matrix dimensions must agree' );
        end
        if ne2 > ne1
            q3  = repmat( quaternion, si2 );
        else
            q3  = repmat( quaternion, si1 );
        end
        for iel = 1 : max( ne1, ne2 )
            q3(iel) = product( q1(min(iel,ne1)), q2(min(iel,ne2)) );
        end
    end % times()
    
    %---------------------------------------------------
    function qm = uminus( q )
        qm  = reshape( quaternion( -double( q )), size( q ));
    end % uminus()
    
%     %---------------------------------------------------
%     function qu = unwrap( q , dim )
%         warning('This may give incorrect results')
%         % Need to find out what the actual process for this is and test it.
%         if ~exist('dim','var')
%             dim = 1;
%         end
%         d = diff(q,1,dim);
%         d = double( q );
%         a = d(2:4,[1 1:end-1],:);
%         b = d(2:4,[1:end],:);
%         norm = @(x)sum(x.^2,1);
%         ia = atan2(norm(cross(a,b)),dot(a,b)); % included angle
%         k = find(ia>pi/2);
%         flip = false(size(ia));
%         for j = 1:numel(k)
%             flip(k(j):end) = ~flip(k(j):end);
%         end
%         sgn = ones(size(flip));
%         sgn(flip) = -1;
%         sgn = shiftdim(sgn,1);
%         qu = q.*sgn;
%     end %unwrap()
    
    %---------------------------------------------------
    %---------------------------------------------------
    
    function [angle, axis] = angleaxis( q )
        % Construct angle-axis pairs equivalent to quaternion rotations
        % Syntax:
        %   [angle, axis] = angleaxis( q )
        %   [angle, axis] = q.angleaxis
        % Input:
        %   q        quaternion array
        % Outputs:
        %   angle    rotation angles in radians
        %   axis     3xN or Nx3 rotation axis unit vectors
        % Note: angle and axis are constructed so at least 2 out of 3
        % elements of axis are >= 0.
        siz         = size( q );
        [angle, s]  = deal( zeros( siz ));
        axis        = zeros( [3 siz] );
        nel         = prod( siz );
        if nel == 0
            return;
        end
        [q, n]      = normalize( q );
        d           = double( q );
        angle(1:end)= 2 * acos( d(1,:) );
        s(1:end)    = sin( 0.5 * angle );
        angle(n==0) = 0;
        s(s==0)     = 1;
        s3          = repmat( shiftdim( s, -1 ), 3, 1 );
        axis(1:end) = reshape( d(2:4,:), [3 siz] )./ s3;
        axis(1,(mod(angle,2*pi)==0)) = 1;
        angle       = chop( angle );
        axis        = chop( axis );
    end %angleaxis()
        
    %---------------------------------------------------
    function v2 = rotate( q, v )
        % Rotate points or vectors
        %   xyz2 = q.rotate(xyz)
        %   xyz2 = rotate(q,xyz)
        %
        % xyz is N-by-3 list of points or vectors which are rotated by q
        % 
        % Rotation is performed using rotation matrices instead of
        % quaternions, as this is up to 7 times faster
        %
        % This function is a simplified version of Mark Tincknell's
        % RotateVector
        
        assert( isa( q, 'quaternion' ), 'q must be a quaternion' )
        assert( numel(q) == 1, 'q must be a single quaternion' )
        assert( size(v,2) == 3, 'v must be N-by-3' )
        R = q.rotationmatrix;
        
        % Convert v from M-by-3-by-N-by... into N-by-3, concatenated 2-D
        % matrix:
        order = 1 : ndims( v ); % [1,2,3,4,...]
        order(end+1) = 2;       % [1,2,3,4,...,2]
        order(2) = [];          % [1,3,4,...,2]     permution order
        vp = reshape( permute(v,order), [], 3 ); 
        
        % Rotate:
        v2 = (R*vp')';
        
        % Reshape back to original form:
        siz = size( v ); 
        order = 1:ndims( v );        % [1,2,3,4,...]
        order([1 2]) = order([2 1]); % [2,1,3,4,...]
        v2 = permute(reshape(v2',siz(order)),order);
    end %Rotate()
    
    %---------------------------------------------------
    function R = rotationmatrix( q )
        % Represent q as a rotation (or direction cosine) matrix
        %
        % Syntax:
        %       R = rotationmatrix(q)
        %       R = q.rotationmatrix
        % Input:
        %       q:  Quaternion array
        % Output:
        %       R:  3-by-3-by-N rotation (or direction cosine) matrix
        %
        % See the following link for solution:
        % http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
        siz = size( q );
        nel = prod( siz );
        R = zeros( [3 3 siz] );
        q = normalize( q );
        for j = 1 : nel
            w = q(j).e(1);
            x = q(j).e(2);
            y = q(j).e(3);
            z = q(j).e(4);
            R(:,:,j) = ...
                [ 1-2*y^2-2*z^2, 2*x*y-2*z*w,   2*x*z+2*y*w
                  2*x*y+2*z*w,   1-2*x^2-2*z^2, 2*y*z-2*x*w
                  2*x*z-2*y*w,   2*y*z+2*x*w,   1-2*x^2-2*y^2 ];
        end
        R = chop( R );
    end %rotationmatrix()
    
end %methods

% Static methods
methods(Static)
    
    
    %---------------------------------------------------
    function q = AngleAxis( angle, axis )
        % Construct quaternions from rotation axes and rotation angles
        % Syntax:
        %   q = quaternion.AngleAxis( angle, axis )
        % Inputs:
        %   angle    array of rotation angles in radians
        %   axis     3xN or Nx3 array of axes (need not be unit vectors)
        % Output:
        %   q        quaternion array
        sig = size( angle );
        six = size( axis );
        [axis, dim, perm] = finddim( axis, 3 );
        if dim == 0
            error( 'axis must have a dimension of size 3' );
        end
        neg = prod( sig );
        nex = prod( six )/ 3;
        if neg == 1
            siz     = six;
            siz(dim)= 1;
            nel     = nex;
            angle   = repmat( angle, siz );
        elseif nex == 1
            siz     = sig;
            nel     = neg;
            axis    = repmat( axis, [1 siz] );
        elseif nex == neg
            siz     = sig;
            nel     = neg;
        else
            error( 'angle and axis must have compatible sizes' );
        end
        for iel = nel : -1 : 1
            d(:,iel) = AngAxis2e( angle(iel), axis(:,iel) );
        end
        q   = quaternion( d );
        q   = reshape( q, siz );
        if neg == 1
            q   = ipermute( q, perm );
        end
    end % quaternion.AngleAxis()
    
    
    %---------------------------------------------------
    function q = EulerAngles( varargin )
        % Construct quaternions from triplets of axes order definitions and
        % Euler angles.
        %
        % Syntax:
        %   q = quaternion.EulerAngles( axes, angles )
        %   q = quaternion.EulerAngles( axes, ang1, ang2, ang3 )
        % Inputs:
        %   axes                string array or cell string array
        %                       '123' = 'xyz' = 'XYZ' = 'ijk', etc.
        %   angles              3xN or Nx3 array of angles in radians  OR
        %   ang1, ang2, ang3    arrays of angles in radians
        % Output:
        %  q                    quaternion array
        
        ics = cellfun( @iscellstr, varargin );
        ic  = cellfun( @ischar, varargin );
        if any( ic )
            varargin{ic} = cellstr( varargin{ic} );
            ics = ic;
        end
        siv     = cellfun( @size, varargin, 'UniformOutput', false );
        axes    = varargin{ics};
        six     = siv{ics};
        nex     = prod( six );

        if nargin == 2  % angles is 3xN or Nx3 array
            angles  = varargin{~ics};
            sig     = siv{~ics};
            [angles, dim, perm] = finddim( angles, 3 );
            if dim == 0
                error( 'Must supply 3 Euler angles' );
            end
            sig(dim)    = 1;
            neg         = prod( sig );
            if nex == 1
                siz     = sig;
                axes    = repmat( axes, siz );
            elseif neg == 1
                siz     = six;
                angles  = repmat( angles, [1 siz] );
            elseif nex == neg
                siz     = sig;
            end
            nel = prod( siz );
            for iel = nel : -1 : 1
                q(iel)  = EulerAng2q( axes{iel}, angles(:,iel) );
            end
    
        elseif nargin == 4  % each of 3 angles is separate input argument
            angles  = conformsize( varargin{~ics} );
            sig     = size( angles{1} );
            neg     = prod( sig );
            if nex == 1
                siz     = sig;
                axes    = repmat( axes, siz );
            elseif neg == 1
                siz     = six;
                for i0 = 1 : 3
                    angles{i0} = repmat( angles{i0}, [1 siz] );
                end
            elseif nex == neg
                siz     = sig;
                axes    = reshape( axes, siz );
            end
            nel = prod( siz );
            for iel = nel : -1 : 1
                q(iel)  = EulerAng2q( axes{iel}, ...
                    [angles{1}(iel), angles{2}(iel), angles{3}(iel)] );
            end
        else
            error( 'Must supply either 2 or 4 input arguments' );
        end % if nargin

        q   = reshape( q, siz );
        if exist( 'perm', 'var' ) && isequal( siz, sig )
            q   = ipermute( q, perm );
        end
        if ~ismatrix( q ) && (size( q, 1 ) == 1)
            q   = shiftdim( q, 1 );
        end
    end % quaternion.EulerAngles()

    %---------------------------------------------------
    function q = nan( varargin )
        % function q = quaternion.nan( siz )
        if isempty( varargin )
            siz = [1 1];
        elseif numel( varargin ) > 1
            siz = [varargin{:}];
        elseif isempty( varargin{1} )
            siz = [0 0];
        elseif numel( varargin{1} ) > 1
            siz = varargin{1};
        else
            siz = [varargin{1} varargin{1}];
        end
        if prod( siz ) == 0
            q   = reshape( quaternion.empty, siz );
        else
            q   = quaternion( nan(siz), nan, nan, nan );
        end
    end % quaternion.nan()

    %---------------------------------------------------
    function q = NaN( varargin )
        % function q = quaternion.NaN( siz )
        q   = quaternion.nan( varargin{:} );
    end % quaternion.NaN()

    %---------------------------------------------------
    function q = ones( varargin )
        % function q = quaternion.ones( siz )
        if isempty( varargin )
            siz = [1 1];
        elseif numel( varargin ) > 1
            siz = [varargin{:}];
        elseif isempty( varargin{1} )
            siz = [0 0];
        elseif numel( varargin{1} ) > 1
            siz = varargin{1};
        else
            siz = [varargin{1} varargin{1}];
        end
        if prod( siz ) == 0
            q   = reshape( quaternion.empty, siz );
        else
            q   = quaternion( ones(siz), 0, 0, 0 );
        end
    end % quaternion.ones()

    %---------------------------------------------------
    function q = rand( varargin )
        % function q = quaternion.rand( siz )
        % Input:
        %  siz      size of output array q
        % Output:
        %  q        uniform random quaternions, normalized to 1,
        %           0 <= q.e(1) <= 1, -1 <= q.e(2:4) <= 1
        if isempty( varargin )
            siz = [1 1];
        elseif numel( varargin ) > 1
            siz = [varargin{:}];
        elseif isempty( varargin{1} )
            siz = [0 0];
        elseif numel( varargin{1} ) > 1
            siz = varargin{1};
        else
            siz = [varargin{1} varargin{1}];
        end
        d   = [ rand( [1, siz] ); 2 * rand( [3, siz] ) - 1 ];
        q   = quaternion( d );
        q   = normalize( q );
        q   = reshape( q, siz );
    end % quaternion.rand()
    
    %---------------------------------------------------
    function q = RotateUtoV( u, v )
        % function q = quaternion.RotateUtoV( u, v ) 
        % Construct quaternions to rotate vectors u into 
        % directions of vectors v 
        % Inputs: 
        %  u, v     3x1 or 3xN or 1x3 or Nx3 arrays of vectors
        % Output:
        %  q        quaternion array
        [u, dimu, permu] = finddim( u, 3 );
        if dimu == 0
            error( 'u must have a dimension of size 3' );
        end
        siu     = size( u );
        siu(1)  = 1;
        neu     = prod( siu );
        [v, dimv, permv] = finddim( v, 3 );
        if dimv == 0
            error( 'v must have a dimension of size 3' );
        end
        siv     = size( v );
        siv(1)  = 1;
        nev     = prod( siv );
        if neu == nev
            siz  = siu;
            nel  = neu;
            perm = permu;
        elseif (neu > 1) && (nev == 1)
            siz  = siu;
            nel  = neu;
            perm = permu;
        elseif (neu == 1) && (nev > 1)
            siz  = siv;
            nel  = nev;
            perm = permv;
        else
            error( 'Number of 3 element vectors in u and v must be 1 or equal' );
        end
        for iel = nel : -1 : 1
            q(iel)  = UV2q( u(:,min(iel,neu)), v(:,min(iel,nev)) );
        end
        q   = ipermute( reshape( q, siz ), perm );
    end % quaternion.RotateUtoV()    
    
    %---------------------------------------------------
    function q = RotationMatrix( R )
        % Construct quaternions from rotation (or direction cosine) matrices
        % Syntax:
        %   q = quaternion.RotationMatrix( R )
        % Input:
        %   R        3x3xN rotation (or direction cosine) matrices
        % Output:
        %   q        quaternion array
        siz = [size(R) 1 1];
        if ~all( siz(1:2) == [3 3] ) || ...
                (abs( det( R(:,:,1) ) - 1 ) > 1e-6 )
            error( 'Rotation matrices must be 3x3xN with det(R) == 1' );
        end
        nel = prod( siz(3:end) );
        for iel = nel : -1 : 1
            d(:,iel) = RotMat2e( chop( R(:,:,iel) ));
        end
        q   = quaternion( d );
        q   = normalize( q );
        q   = reshape( q, siz(3:end) );
    end % quaternion.RotationMatrix()
    
    %---------------------------------------------------
    function q = zeros( varargin )
        % Construct zero quaternions
        %   q = quaternion.zeros( siz )
        if isempty( varargin ) || isempty( varargin{1} )
            siz = [1 1];
        elseif numel( varargin ) > 1
            siz = [varargin{:}];
        elseif numel( varargin{1} ) > 1
            siz = varargin{1};
        else
            siz = [varargin{1} varargin{1}];
        end
        q   = quaternion( zeros(siz), 0, 0, 0 );
    end % quaternion.zeros()
    
end %methods(Static)

end %classdef quaternion
                
  
                
                
                
% ========================================================================
%   Helper functions

% ------------------------------------------------------------------------
function eout = AngAxis2e( angle, axis )
% Convert one Angle-Axis -> one quaternion
% Syntax:
%   eout = AngAxis2e( angle, axis )
s   = sin( 0.5 * angle );
v   = axis(:);
vn  = norm( v );
if vn == 0
    if s == 0
        c   = 0;
    else
        c   = 1;
    end
    u   = zeros( 3, 1 );
else
    c   = cos( 0.5 * angle );
    u   = v(:) ./ vn;
end
eout    = chop( [ c; s * u ] );
if (eout(1) < 0) && (mod( angle/(2*pi), 2 ) ~= 1)
    eout = -eout; % rotationally equivalent quaternion with real element >= 0
end
end %AngAxis2e()


% ------------------------------------------------------------------------
function qout = EulerAng2q( axes, angles )
% function qout = EulerAng2q( axes, angles )
% One triplet Euler Angles -> one quaternion
na   = length( axes );
axis = zeros( 3, na );
for i0 = 1 : na
    switch axes(i0)
        case {'1', 'i', 'x', 'X'}
            axis(:,i0) = [ 1; 0; 0 ];
        case {'2', 'j', 'y', 'Y'}
            axis(:,i0) = [ 0; 1; 0 ];
        case {'3', 'k', 'z', 'Z'}
            axis(:,i0) = [ 0; 0; 1 ];
        otherwise
            error( 'Illegal axis designation' );
    end
end
q0  = quaternion.AngleAxis( angles(:).', axis );
for i0 = 1 : numel( q0 )
    if i0 == 1
        qout = q0(i0);
    else
        qout = q0(i0) * qout;
    end
end
if qout.e(1) < 0
    qout = -qout; % rotationally equivalent quaternion with real element >= 0
end
end % EulerAng2q()


% ------------------------------------------------------------------------
function eout = RotMat2e( R )
% function eout = RotMat2e( R )
% One Rotation Matrix -> one quaternion
if all( all( R == 0 ))
    eout    = zeros(4,1);
else
    eout(1) = 0.5 * sqrt( max( 0, R(1,1) + R(2,2) + R(3,3) + 1 ));
    if eout(1) == 0
        eout(2) = sqrt( max( 0, -0.5 *( R(2,2) + R(3,3) ))) * ...
                        sgn( -R(2,3) );
        eout(3) = sqrt( max( 0, -0.5 *( R(1,1) + R(3,3) ))) * ...
                        sgn( -R(1,3) );
        eout(4) = sqrt( max( 0, -0.5 *( R(1,1) + R(2,2) ))) * ...
                        sgn( -R(1,2) );
    else
        eout(2) = 0.25 *( R(3,2) - R(2,3) )/ eout(1);
        eout(3) = 0.25 *( R(1,3) - R(3,1) )/ eout(1);
        eout(4) = 0.25 *( R(2,1) - R(1,2) )/ eout(1);
    end
end
eout    = chop( eout(:) );
if eout(1) < 0
    eout = -eout; % rotationally equivalent quaternion with real element >= 0
end
end % RotMat2e()


% ------------------------------------------------------------------------
function s = sgn( x )
% function s = sgn( x ), if x >= 0, s = 1, else s = -1
s   = ones( size( x ));
s(x < 0) = -1;
end % sgn


% ------------------------------------------------------------------------
function qout = UV2q( u, v )
% function qout = UV2q( u, v )
% One pair vectors U, V -> one quaternion
w       = cross( u, v );    % construct vector w perpendicular to u and v
magw    = norm( w );
dotuv   = dot( u, v );
if magw == 0
% Either norm(u) == 0 or norm(v) == 0 or dotuv/(norm(u)*norm(v)) == 1
    if dotuv >= 0
        qout    = quaternion( [ 1; 0; 0; 0 ] );
        return;
    end
% dotuv/(norm(u)*norm(v)) == -1
    magv    = norm( v );
% If v == [v(1); 0; 0], rotate by pi about the [0; 0; 1] axis
    if magv == abs( v(1) )
        qout    = quaternion( [ 0; 0; 0; 1 ] );
        return;
    end
% Otherwise constuct "what" such that dot(v,what) == 0, and rotate about it
% by pi
    w2      = v(3) / sqrt( v(2)^2 + v(3)^2 ); 
    what    = [ 0; w2; sqrt( 1 - w2^2 ) ];
    costh   = -1;
else
% Use w as rotation axis, angle between u and v as rotation angle
    what    = w(:) / magw;
    costh   = dotuv /( norm(u) * norm(v) );
end
c       = sqrt( 0.5 *( 1 + costh ));    % real element >= 0
s       = sqrt( 0.5 *( 1 - costh ));
eout    = [ c; s * what ];
qout    = quaternion( chop( eout ));
end % UV2q()


% ------------------------------------------------------------------------
function varargout = qcheck(opt,varargin)
% Check input(s) are quaternions, if not either error or assign
isquat = cellfun(@(x)isa(x,'quaternion'),varargin);
varargout = varargin;   % Initialize
if all(isquat)
    return
end
switch opt
    case 'error'
        msg = 'Input(s) must be quaternions';
        ds  = dbstack;
        msgid = [ds(2).name ':NotQuaternion'];
        throwAsCaller(MException(msgid,msg))
    case 'assign'
        try
        idx = find(~isquat);
        for j = idx
            vj = varargin{j};
            varargout{j} = quaternion(real(vj),imag(vj));
        end
        catch ME
            throwAsCaller(ME)
        end
end %switch
end % qcheck()


% ------------------------------------------------------------------------
function out = chop( in, tol )
% Replace values that differ from an integer by <= tol by the integer
% Inputs:
%   in       input array
%   tol      tolerance, default = eps
% Output:
%   out      input array with integer replacements, if any
if ~exist( 'tol', 'var' ) || isempty( tol )
    tol = eps;
end
out = in;
rin = round( in );
lx  = abs( rin - in ) <= tol;
out(lx) = rin(lx);
end % chop()


% ------------------------------------------------------------------------
function out = conformsize( varargin )
% Replicate scalar arguments to the size of argument with largest number of
% elements, and reshape arrays with the same number of elements
% Inputs:
%  arg1, arg2, ...  arrays of the same size or scalars
% Output:
%  out              cell array of input arrays or expanded scalars
nelem     = cellfun( 'prodofsize', varargin );
[nel, nei] = max( nelem );
siz       = size( varargin{nei} );
for i0 = nargin : -1 : 1
    if nelem(i0) == 1
        out{i0} = repmat( varargin{i0}, siz );
    elseif nelem(i0) == nel
        out{i0} = reshape( varargin{i0}(:), siz );
    elseif ~isequal( siz, size( varargin{i0} ))
        error( 'Inputs must be the same sizes or scalars' );
    else
        out(i0) = varargin(i0);
    end
end
end % conformsize()


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
    perm = [ dim : ndm, 1 : dim-1 ];
    aout = permute( ain, perm );
end
end % finddim()


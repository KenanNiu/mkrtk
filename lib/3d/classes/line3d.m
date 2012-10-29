classdef line3d
% LINE3D class.
%
%
% Constructors:
%   l = line3d( point3d1, point3d2 )
%   l = line3d( point3d,  vector3d )
%   l = line3d( 'point-point', point1, point2 )
%   l = line3d( 'point-vector', point, direction )
%
% Methods:
%
%
% Accessing Properties:
%
%
% Internal Implementation:
%
%
%
% See also POINT3D, VECTOR3D

properties (SetAccess = protected, GetAccess = protected)
    p
    v
    p1
    p2
end


methods % ordinary methods
    function l = line3d( varargin ) %(constructor)
        switch nargin
            case 0 % nargin == 0; called by higher nargin cases
                pe = point3d.empty;
                l.p  = pe;
                l.v  = vector3d.empty;
                l.p1 = pe;
                l.p2 = pe;
                return
                
            case 2 % line from 2 objects:
                [arg1,arg2] = deal(varargin{:});
                if isa(arg1, 'point3d') && isa(arg2, 'vector3d')
                    arg2 = arg2.normalize;
                    siz = size(arg1);
                    for j = numel(arg1) : -1 : 1
                        l(j).p  = arg1(j);
                        l(j).v  = arg2(j);
                    end
                elseif isa(arg1, 'point3d') && isa(arg2, 'point3d')
                    siz = size(arg1);
                    for j = numel(arg1) : -1 : 1
                        l(j).p1 = arg1(j);
                        l(j).p2 = arg2(j);
                    end
                else
                    error(['Invalid inputs for constructing with two argments. ',...
                        'You provided a %s and a %s, where inputs should be ',...
                        'a POINT3D and a VECTOR3D, or 2x POINT3D objects'], ...
                        upper(class(arg1)), upper(class(arg2)))
                end
                
            case 3 % line from doubles
                [opt,arg1,arg2] = deal(varargin{:});
                switch lower(opt)
                    case 'point-vector'
                        pt = point3d(arg1);
                        vec = vector3d(arg2).normalize;
                        siz = size(pt);
                        for j = numel(pt) : -1 : 1
                            l(j).p = pt(j);
                            l(j).v = vec(j);
                        end
                    case 'point-point'
                        pt1 = point3d(arg1);
                        pt2 = point3d(arg2);
                        siz = size(pt1);
                        for j = numel(pt1) : -1 : 1
                            l(j).p1 = pt1(j);
                            l(j).p2 = pt2(j);
                        end
                    otherwise
                        error(['Invalid constructor option %s. Use either ',...
                            'POINT-POINT, or POINT-VECTOR.'], upper(opt))
                end %switch
                
            otherwise
                error('Incorrect number of input arguments')
        end %switch
        l = reshape(l, siz );
    end %line3d()
    
    %----------------------------------------------------
    function v = direction(l)
        %DIRECTION Generalised method for getting the 
        %  direction vector of a line, regardless of what 
        %  mode it is defined in.
        v = l.toray.vector; 
    end %direction()
    
    %----------------------------------------------------
    function display(l)
        %DISPLAY Called during unsupressed output to the
        % command line.
        % See also disp() for standard un-formatted output.
        siz = size( l );
        nel = [1 cumprod( siz )];
        ndim = length(siz);
        disp(' ')
        fprintf('%s = \n\n',inputname(1));
        sstr = sprintf('%dx',siz);  sstr(end) = [];
        fprintf('  %s <a href="matlab:help line3d">line3d</a>',sstr);
        fprintf(' with <a href="matlab:methods(''line3d'')">methods</a>\n\n');
        for iel = 1 : nel(end)
            if nel(end) == 1
                sub = '';
            else
                sub = ')\t';
                jel = iel - 1;
                for idm = ndim : -1 : 1
                    idx = floor( jel / nel(idm) ) + 1;
                    sub = [',' int2str(idx) sub]; %#ok<AGROW>
                    jel = rem( jel, nel(idm) );
                end
                sub(1) = '(';
            end
            if isempty(l(iel).p2)
                arg1name = 'point';
                arg2name = 'vector';
                arg1val  = [l(iel).point.xyz];  % brackets handle empty results
                arg2val  = [l(iel).vector.uvw];
                arg1cls  = class(l(iel).point);
                arg2cls  = class(l(iel).vector);
            else
                arg1name = 'point1';
                arg2name = 'point2';
                arg1val  = [l(iel).point1.xyz];
                arg2val  = [l(iel).point2.xyz];
                arg1cls  = class(l(iel).point1);
                arg2cls = class(l(iel).point2);
            end
            if isempty(arg1val)
                arg1val = sprintf('[empty %s]',arg1cls);
            else
                arg1val = sprintf('[%g %g %g]',arg1val);
            end
            if isempty(arg2val)
                arg2val = sprintf('[empty %s]',arg2cls);
            else
                arg2val = sprintf('[%g %g %g]',arg2val);
            end
            fprintf( ['%s' sub ': %s: %s\t%s: %s\n'], ...
                inputname(1), arg1name, arg1val, arg2name, arg2val )
        end %for
        fprintf('\n');
    end %display()
    
    %----------------------------------------------------
    function tf = eq(l1,l2)
        %EQ Test equality of elements.
        % Usage:
        %       l1 == l2
        %       l1.eq(l2)
        %
        % This function can be used when l1 and/or l2 contain mixed modes
        % (ie, some line3d objects in the array are point-vector, and some
        % are point-point). 
        %
        % Either/both l1 and l2 can be arrays or 1-by-1.  In the case where
        % a scalar is compared with an array, it would be more consise to
        % use repmat() then compare two arrays of equal size, but it is
        % faster to handle them explicitly, which is what the nested IF
        % statement does.
        if strcmpi( class(l1), class(l2) ) 
            si1 = size( l1 );
            si2 = size( l2 );
            ne1 = prod( si1 );
            ne2 = prod( si2 );
            l1isseg = l1.issegment;
            l2isseg = l2.issegment;
            if (ne1 == 0) || (ne2 == 0)
                % Both empty
                tf   = logical([]);
                return;
            elseif ne1 == 1
                % Compare 1-by-1 with N-by-M-by...
                tf = false(size(l2));
                modeeq = l1isseg == l2isseg;
                switch l1isseg
                    case 0
                        tf(modeeq) = ( l1.p == [l2(modeeq).p] ) & ( l1.v == [l2(modeeq).v] );
                    case 1
                        tf(modeeq) = ( l1.p1 == [l2(modeeq).p1] ) & ( l1.p2 == [l2(modeeq).p2] );
                end
            elseif ne2 == 1
                % Compare N-by-M-by... with 1-by-1
                tf = false(size(l1));
                modeeq = l1isseg == l2isseg;
                switch l2isseg
                    case 0
                        tf(modeeq) = ([l1(modeeq).p] == l2.p) & ([l1(modeeq).v] == l2.v);
                    case 1
                        tf(modeeq) = ([l1(modeeq).p1] == l2.p1) & ([l1(modeeq).p2] == l2.p2);
                end
            elseif isequal( si1, si2 )
                % Compare two equally sized arrays
                tf = false(size(l1));
                rays = ~l1isseg & ~l2isseg;
                segs = l1isseg & l2isseg;
                tf(rays) = ([l1(rays).p]  == [l2(rays).p] ) & ([l1(rays).v]  == [l2(rays).v] );
                tf(segs) = ([l1(segs).p1] == [l2(segs).p1]) & ([l1(segs).p2] == [l2(segs).p2]);
            else
                error('Matrix dimensions must agree')
            end
        else
            tf = false;
        end        
    end %eq()
    
    %----------------------------------------------------
    function o = origin(l)
        % ORIGIN Collect origin of line
        %   Rays     => p  (.point)
        %   Segments => p1 (.point1)
        o = repmat(point3d(),size(l));  % Initialize
        isseg = l.issegment;
        isray = ~isseg;
        if any(isseg(:))
            o(isseg) = l(isseg).point1;     
        end
        if any(isray(:))
            o(isray) = l(isray).point;
        end     
    end %origin()
    
    %----------------------------------------------------
    function tf = ne(l1,l2)
        tf = ~eq( l1, l2 );
    end %ne()
    
    %----------------------------------------------------
    %function tf = isfinite(l)
    %end
    
    %----------------------------------------------------
    function tf = isray(l)
        % ISRAY Test if this line is an infinite line, or ray
        %   (ie, not a finite segment)
        %--- comprehensive version:
        %ep  = cellfun(@isempty, {l.p} );
        %ev  = cellfun(@isempty, {l.v} );
        %ep1 = cellfun(@isempty, {l.p1});
        %ep2 = cellfun(@isempty, {l.p2});
        %tf = ~ep & ~ev & ep1 & ep2;         % p & v not empty; p1 & v1 empty
        %--- fast version:
        tf = ~cellfun(@isempty, {l.v});
        %---
        tf = reshape(tf,size(l));
    end %isray()
    
    %----------------------------------------------------
    function tf = issegment(l)
        % ISSEGMENT Test if this line is a finite line segment
        %   (ie, not an infinite line or ray)
        tf = ~isray(l);
    end %issegment()
    
    %----------------------------------------------------
    function i = project(l,obj)
        % PROJECT Project an object (point3d only?) onto line
        if numel(l) > numel(obj)
            siz = size(l);
        else
            siz = size(obj);
        end
        switch lower( class(obj) )
            case 'point3d'
                l = l.toray;
                dp = bsxfun(@minus, [obj.xyz]', [l.point.xyz]');
                dv = [l.vector.uvw]';
                i  = bsxfun(@rdivide, sum(bsxfun(@times, dp, dv), 2), sum(dv.^2, 2));
                % i is a scalar
            otherwise
                error('Unhandled object')
        end
        i = reshape(i,siz);
    end %project()
    
    %----------------------------------------------------
    function lout = toray(l)
        % TORAY Convert any point-point lines in LIN to 
        %   point-vector
        lout = l;
        idx = find(l.issegment);
        if ~isempty(idx)
            pt1 = l(idx).point1;
            pt2 = l(idx).point2;
            vec = vector3d([pt2.xyz] - [pt1.xyz]).normalize;
            lout(idx) = line3d(pt1,vec);
        end
    end %toray()
    

    % ===================================================
    % The following methods act similarly to getter methods.
    % However, we do not use getter methods because they act 
    % on each individual object, rather than on an array of
    % objects.  We want to act on the array of objects and
    % return arrays as outputs.
    %
    % Thus in these functions need to adopt names different
    % to the property names, and we also protect the 
    % properties.
    % Property / get method correspondence is:
    %      p:  point()
    %      v:  vector()
    %     p1:  point1()
    %     p2:  point2()
    %
    
    
    %----------------------------------------------------
    function pt = point(l)
        pt = [l.p];
        if (numel(pt) ~= numel(l))
            assert(numel(pt) == numel(l), get_prop_errstr(),'point')
        else
            pt = reshape(pt, size(l));
        end
    end %point()
    
    %----------------------------------------------------
    function vec = vector(l)
        vec = [l.v];
        if (numel(vec) ~= numel(l))
            assert(numel(vec) == numel(l), get_prop_errstr(),'vector')
        else
            vec = reshape(vec, size(l));
        end
    end %vector()
    
    %----------------------------------------------------
    function pt1 = point1(l)
        pt1 = [l.p1];
        if (numel(pt1) ~= numel(l))
            assert(numel(pt1) == numel(l), get_prop_errstr(),'point1')
        else
            pt1 = reshape(pt1, size(l));
        end
    end %point1()
    
    %----------------------------------------------------
    function pt2 = point2(l)
        pt2 = [l.p2];
        if (numel(pt2) ~= numel(l))
            assert(numel(pt2) == numel(l), get_prop_errstr(),'point2')
        else
            pt2 = reshape(pt2, size(l));
        end
    end %point2()
    
end %methods
    
end %classdef

% ========================================================================
%   Helper functions


% ------------------------------------------------------------------------
function msg = get_prop_errstr
msg = ['You cannot mass-collect the %s of a set of line3d ',...
    'objects that are defined in different modes (point-vector & ',...
    'point-point).'];
end %get_propo_ermsg()


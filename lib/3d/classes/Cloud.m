classdef Cloud
% Cloud class

properties
    Tag = ''
    Path = ''
    xyz
    Origin
    CS
    Visible = true
end

properties (Hidden=true)
    traceLengths = [];
end

properties (Hidden=true, SetAccess=private)
    Version = 2;        
end

methods
    % ------------------------------------------
    function c = Cloud( varargin )
        % Cloud()       1-by-1 empty Cloud
        % Cloud([])     0-by-0 Cloud
        % Cloud(xyz)    1-by-1 Cloud with points, origin, and CS configured
        % Cloud(X)    
        % Cloud(struct) 
        % Cloud(...,tag)
        if nargin == 0
            % Cloud()
            return
        end
        if isempty( varargin{1} )
            c = Cloud.empty;
            return
        end
        
        if isa(varargin{1},'struct')
            c = structs2clouds( c, varargin{:} );
        else
            c = data2cloud( c, varargin{:} );
        end
        

    end %Cloud()
    
    % ------------------------------------------
    function c = downsample(c,v)
        % c.downsample( int32(150) )    % downsample to N points
        % c.downsample( 0.5 )           % downsample to min spacing
        % 
        % Worth considering an option for display or rendering using the
        % Ramer Douglas Peucker algorithm:
        %   http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm 
        
        if numel(v) == 1
            if isinteger( v )
                % Simplest, most rudimentary downsample to n points
                nnew = v;
                np = size(c.xyz,1);
                if nnew >= np
                    return
                end
                idx = round( linspace(1,np,nnew) );
                sel = false(np,1);
                sel(idx) = true;
            else
                % downsample to minimum spacing
                %ds = v;
                error('not yet implemented')
            end
        elseif islogical(v)
            assert(numel(v) == size(c.xyz,1),...
                'Logical index vector is incorrectly sized.')
            sel = v;
        else
            error('Unhandled input')                
        end
        
        % Downsample:
        c.xyz = c.xyz(sel,:);
        
        % Adjust trace partitioning if necessary:
        if ~isempty(c.traceLengths)
            tfc = mat2cell(sel,c.traceLengths);
            c.traceLengths = cellfun(@sum,tfc);
        end
    end %downsample()
    
    % ------------------------------------------
    function [r,t] = gettransform(c1,c2,rdesc)
        % Get transformation that maps C1->C2 as described in the following
        % pseudo code:
        %   C2 = R*C1 + T
        %
        % Usage:
        %   [R,T] = gettransform(c1,c2,'rotmat')  
        %   [Q,T] = gettransform(c1,c2,'quat')
        %
        % Ouputs: If C1 & C2 contain N elements, the outputs are sized as
        % follows:
        %   R => 3-by-3-by-N
        %   Q => N-by-1
        %   T => N-by-3
        assert( numel(c1)==numel(c2), ...
            'C1 and C2 must have the same number of elements')
        
        % Calculate transformation(s):
        for j = numel(c1):-1:1
            r(:,:,j) = c2(j).CS*c1(j).CS';                              % R = CS2'*CS1
            t(j,:) = ( c2(j).Origin(:) - r(:,:,j)*c1(j).Origin(:) )';   % t = O2 - R*O1
        end
        
        % Manage output form:
        rdesc = lower(rdesc);
        descHas = @(s)~isempty(strfind(rdesc,s));
        if descHas('mat')           % Rotation matrix, fine
            return
        elseif descHas('quat')      % Quaternion, convert
            r = quaternion.RotationMatrix(r);
            r = reshape(r,size(c1));
        else                        % Other, error
            error('Unknown rotation description: %s',upper(rdesc))
        end
    end %gettransform()
    
    % ------------------------------------------
    function varargout = plot(c,varargin)
        % c.plot
        % c.plot(axs)
        % c.plot(...,Property,Value,Property,Value,...)
        [hax,propVal] = axisHandleFromVarargin(varargin);
        h = zeros(size(c));
        for j = 1:numel(c);
            h(j) = plot3(c(j).xyz(:,1),c(j).xyz(:,2),c(j).xyz(:,3),...
                propVal{:},...
                'Parent',hax,...
                'DisplayName',c(j).Tag,...
                'Tag',c(j).Tag);
        end
        if nargout > 0
            varargout = {h};
        end
    end %plot()
    
    % ------------------------------------------
    function varargout = plotcs(c,varargin)
        % c.plotcs
        % c.plotcs(hax);
        hax = axisHandleFromVarargin(varargin);
        hcs = zeros(size(c));
        for j = 1:numel(c)
            scale = max( std(c(j).xyz) );
            hcs(j) = drawCS(hax,c(j).Origin,c(j).CS,scale,c(j).Tag);
        end
        if nargout > 0
            varargout = {hcs};
        end
    end %plotcs()
    
    % ------------------------------------------
    function varargout = plottraces(c,varargin)
        % c.plottraces
        % c.plottraces(ax)
        % c.plottraces(...,Property,Value,Property,Value,...)
        % hg = c.plottraces(...)
        [hax,propVal] = axisHandleFromVarargin(varargin);
        ok = ~cellfun(@isempty,{c.traceLengths});
        assert(all(ok),'Cannot plot traces: No traces defined for one or more clouds')
        HOLD = ishold(hax);
        nc = numel(c);
        hg = zeros(nc,1);
        if ~HOLD
            % We turn holding on to plot all traces, then turn it off (if
            % required)
            cla(hax,'reset');
            hold(hax,'on')
        end
        for j = 1:nc
            hg(j) = hggroup('Parent',hax,'Tag',c(j).Tag);
            traces = c(j).traces;
            nt = numel(traces);
            ht = zeros(nt,1);
            for k = 1:numel(traces)
                tk = traces{k};
                ht(k) = plot3(tk(:,1),tk(:,2),tk(:,3),...
                    propVal{:},'Parent',hax,'Parent',hg(j));
            end
        end
        if ~HOLD
            view(3)
            hold(hax,'off')
        end
        if nargout > 0
            varargout = {hg};
        end
    end %plottraces()
    
    % ------------------------------------------
    function c = setorigin(c)
        for j = 1:numel(c)
            c(j).Origin = mean(c(j).xyz);
        end
    end %setorigin()
    
    % ------------------------------------------
    function c = setprincompcs(c)
        for j = 1:numel(c)
            pc = princomp(c(j).xyz);               % Principal components
            pc(:,3) = cross(pc(:,1),pc(:,2));   % Ensure RH coordinate system
            c(j).CS = pc;
        end
    end %setprincompcs()
    
    % ------------------------------------------
    function c = transform(c,r,v)
        % c.transform(R,v)  Performs rotation around the origin by the 
        %                   rotation matrix R, then translation by vector v
        % c.transform(q,v)  Performs rotation around the origin by 
        %                   quaternion q, then translation by vector v
        
        % Handle formats of v:
        if numel(v) == 3
            v = v(:)';
        else
            assert(size(v,2) == 3, ...
                'Translation vector must be 1-by-3 or N-by-3')
        end
        nv = size(v,1);
        % Handle formats of r:
        if isa(r,'quaternion')
            nr = numel(r);
            assert(nr == nv,[...
                'Translation vector V must contain as many rows',...
                ' as quaternion Q'])
        elseif isnumeric(r) && isequal( [3 3], size(r(:,:,1)))
            nr = size(r,3);
            assert(nr == nv, [...
                'Rotation matrix R must contain as many layers as ',...
                'translation vector V contains rows'])
        else
            error(['Bad rotation definition. Define a rotation with ',...
                'either a rotation matrix or a quaternion']);
        end
        % Now rotate + translate:
        n = numel(c);
        if (n == 1) || (n > 1 && nr == 1)
            % Less overhead this way if all transformations are the same:
            c = c.rotate(r).translate(v);
        else
            % Otherwise individually call rotation/translate:
            for j = 1:n
                c(j) = c(j).rotate(r(:,:,j)).translate(v(j,:));
            end
        end
    end %transform()
    
    % ------------------------------------------
    function c = translate(c,t)
        % Translate clouds by specified vector T
        % Usage:
        %   c.translate([x y z])    % Single cloud, single vector   => 1 cloud
        %   C.translate([x y z])    % N clouds, single vector       => N clouds
        %   C.translate(XYZ)        % N clouds, N vectors           => N clouds
        %   c.translate(XYZ)        % single cloud, N vectors       => N clouds
        nc = numel(c);
        nt = size(t,1);
        assert( (ndims(t)==2) && (size(t,2)==3),...
            'Translation vector must be N-by-3') %#ok<ISMAT>
        assert( ...
            (nc == nt) || ...   % includes simple case: nc==1 && nt==1
            (nc == 1)  || ...   % 1 cloud, N translations
            (nt == 1) ,...      % N clouds, 1 translation
            'Input arguments incorrectly sized. See help for details')
        % Now we expand c or t to the correct size
        if nc > 1 && nt == 1
            t = t(ones(nc,1),:);    % replicate (faster than repmat)
        end
        if nt > 1 && nc == 1
            c = c(ones(nt,1),:);
        end
        % Then do element-by-element translations
        for j = 1:numel(c)
            c(j).xyz = shift(c(j).xyz,t(j,:));           % Translate points
            c(j).Origin = shift(c(j).Origin,t(j,:));     % Translate origin
        end
    end %translate()
    
    % ------------------------------------------
    function r = rois(c)
        % Alias for TRACES()
        r = c.traces;
    end %rois
    
    % ------------------------------------------
    function c = rotate(c,varargin)
        % Rotate clouds using a quaternion, or an angle & axis definition
        %
        % Usage - quaternion rotation:
        %   c.rotate(Q)
        %   c.rotate(Q,ORIGIN)
        %
        % Usage - angle & axis:
        %   c.rotate(ALPHA, [THETA PHI])
        %   c.rotate(ALPHA, [X Y Z])
        %   c.rotate(ALPHA, ..., ORIGIN)
        %
        % Usage - rotation matrix:
        %   c.rotate(R)
        %   c.rotate(R,ORIGIN)
        %
        % The format for specifying angle & axis form of rotation is
        % identical to that of the matlab function ROTATE. See ROTATE for
        % further details. 
        sizc = size(c);
        nc   = numel(c);
        rctr = [];
        nvargs = numel(varargin);
        if isa(varargin{1},'quaternion')
            %c.rotate(Q)
            %c.rotate(Q,ORIGIN)
            q = varargin{1};
            if nvargs >= 2
                rctr = varargin{2};
            end
            % Expand clouds or q if required:
            if nc==1 && numel(q)>1
                sizc = size(q);     % Update size
                c = c(ones(sizc));  % Expand clouds
            elseif nc>1 && numel(q)==1
                q = q(ones(sizc));  % Expand quats
            end
            % Rotation function
            rfun = @(i,pts)q(i).rotate(pts);
            
        elseif length(varargin{1}) == numel(varargin{1})    % Angle is 1xN or Nx1
            %c.rotate(ALPHA, [THETA PHI])
            %c.rotate(ALPHA, [THETA PHI], ORIGIN)
            %c.rotate(ALPHA, [X Y Z])
            %c.rotate(ALPHA, [X Y Z], ORIGIN)
            [alfa,azel] = varargin{1:2};
            assert(ismatrix(alfa) && size(alfa,2)==1,...
                'Angle should be N-by-1')
            assert(ismatrix(azel) && (size(azel,2)==2 || size(azel,2)==3),...
                'Rotation angles or axis should be N-by-2 or N-by-3 respectively')
            if nvargs >= 3
                rctr = varargin{3};
            end
            % Expand clouds or alfa & azel if required:
            if nc==1 && numel(alfa)>1
                sizc = size(alfa);  % Update size
                c = c(ones(sizc));  % Expand clouds
            elseif nc>1 && numel(alfa)==1
                rep = ones(sizc);
                alfa = alfa(rep);   %\ Expand rotation
                azel = azel(rep,:); %/             
            end
            % Rotation function:
            rfun = @(i,pts)rotatePts(pts,alfa(i),azel(i,:));
            
        elseif isequal( [3 3], [size(varargin{1},1), size(varargin{1},1)]) % 3x3x...
            %c.rotate(R)
            %c.rotate(R,ORIGIN)
            R = varargin{1};
            nr = size(R,3);
            if nvargs >= 2
                rctr = varargin{2};
            end
            % Expand clouds or R if required
            if nc==1 && nr>1
                sizc = [nr,1];              % Update size
                c = c(ones(sizc));          % Expand clouds
            elseif nc>1 && nr==1
                R = R(:,:,ones(1,1,nc));    % Expand rotation
            end
            % Rotation function
            rfun = @(i,pts)(R(:,:,i)*(pts'))';
            
        else
            error(['Bad rotation definition.  First input should be a ',...
                'quaternion, a rotation angle (followed by a rotation ',...
                'axis), or a rotation matrix.'])
        end
        
        
        % Translate:
        if ~isempty(rctr)
            c = c.translate(-rctr);
        end
        
        % Rotate:
        for j = 1:numel(c)
            c(j).xyz = rfun(j,c(j).xyz);
            c(j).CS  =(rfun(j,c(j).CS'))';
            c(j).Origin = rfun(j,c(j).Origin);
        end
        
        % Translate back
        if ~isempty(rctr)
            c = c.translate(rctr);
        end
        
    end %rotate()
        
    % ------------------------------------------
    function t = traces(c)
        assert(numel(c)==1,'This function cannot be called on an array')
        if isempty(c.traceLengths)
            t = [];
            %error('Cloud was not created from discrete traces')
        else
            t = mat2cell(c.xyz,c.traceLengths);
        end
    end %traces()
        
end %methods


methods(Static)
    
    % ------------------------------------------
    function [clouds,phases] = Load(filename,varargin)
        % Create new cloud(s) frome data contained in file.
        %   File could be either
        %       - Ascii text file containing an N-by-3 matrix of data that
        %         can be read with matlab's LOAD() function
        %       - *.mat file containing ROI objects
        %       - *.mat file containing ROI structures (old versions)
        %
        % c = Cloud.Load(filename)
        % c = Cloud.Load(filename,tag)
        % c = Cloud.Load(filename,tag,N)
        phases = [];
        [~,fname,ext] = fileparts(filename);
        data = load(filename);
        
        % Choose loading method based on file extension & content:
        ISMAT = isequal(ext,'.mat');
        ISASCII = ~ISMAT;
        VERLESSTHAN2 = isfield(data,'xyz');
        ISROIOBJ = isfield(data,'ROI');
        
        if ISASCII
            clouds = Cloud(data,fname);  % since data is N-by-3
        elseif VERLESSTHAN2
            X = data.xyz;   %N-by-1 cell array of traces in 3d
            clouds = Cloud(X);
        elseif ISROIOBJ
            % Standard load: ROI objects with version >= 2
            % In this case we can have multiple phases in the one file, so
            % we need to break them out into separate structures
            [clouds,phases] = roi2cloud(data.ROI);
        else
            error('Unhandled case')
        end
        % Set path:
        [clouds.Path] = deal(filename);
        
        % Set Tag:
        if numel(varargin) == 0  || isempty(varargin{1})     % Set default
            [clouds.Tag] = deal([fname,ext]);
        else                                                 % Set specified
            [clouds.Tag] = deal(varargin{1});       
        end
        
        % Downsample:
        if numel(varargin) >= 2
            clouds = clouds.downsample(N);
        end
        
    end %Load()
    
end %methods(Static)

end %classdef


% ------------------------------------------------------------------------
function [hax,v] = axisHandleFromVarargin(v)
% Get axis handle from varargin, if it exists
% [hax,v] = axisHandleFromVarargin(varargin)
if numel(v) >= 1 && isnumeric( v{1} ) 
    if ishghandle( v{1} ) && strcmpi(get(v{1},'Type'),'axes')
        hax = v{1};
    else
        error('Axis handle is invalid')
    end
    v(1) = [];
else
    hax = gca;
end
end %axisHandleFromVarargin()


% ------------------------------------------------------------------------
function c = data2cloud(c,varargin)
% Populate/expand c with data from cell array or matrix
if isa( varargin{1}, 'cell' )
    X = varargin{1};
    c.traceLengths = cellfun(@(x)size(x,1),X);
    pts = cat(1,X{:});
else
    pts = varargin{1};
end
if size(pts,2)~=3
    pts = pts';
end
if ndims(pts) > 2 || size(pts,2) ~= 3 %#ok<ISMAT>
    error(['First input (cloud points) should be a 2D',...
        ' N-by-3 or 3-by-N matrix, or cell array of',...
        ' N-by-3 matrices'])
end
% Set properties:
c.xyz = pts;
c = c.setorigin;
c = c.setprincompcs;

if numel(varargin) >= 2
    c.Tag = varargin{2};
end
end %data2cloud()


% ------------------------------------------------------------------------
function hcs = drawCS(hax,origin,CS,scale,parent_tag)
%DRAWCS Draw the coordinate system as a hggroup
%
% Inputs:
%          hax: Matlab axes to plot on
%       origin: 1-by-3 Origin
%           CS: 3-by-3 Coordinate system system
%        scale: Scalar, size of the lines plotted
%   parent_tag: Tag of the associated point cloud 
%
% Build a hggroup with the following structure:
%   h_coord_system              (hggroup)
%       |--> vectors            (hggroup)
%       |       |--> line_x        (line)
%       |       |--> line_y        (line)
%       |       |--> line_z        (line)
%       |--> labels             (hggroup)
%               |--> text_x        (text)
%               |--> text_y        (text)
%               |--> text_z        (text)

origin = origin(:)';    % Just check

hold(hax,'on')

hcs = hggroup('Parent',hax,'Tag', [parent_tag ' corrdinate system'],'HitTest','off');
hv  = hggroup('Parent',hcs,'Tag','vectors','HitTest','off');
hl  = hggroup('Parent',hcs,'Tag','labels','HitTest','off');
clr = [0.1, 0, 0.3];

% Axis vectors:
vprops = {'Color',clr,'Parent',hv,'HitTest','off'};
o = {origin(1) origin(2) origin(3)};
qx = CS(:,1)*scale;
qy = CS(:,2)*scale;
qz = CS(:,3)*scale;
quiver3(o{:},qx(1),qx(2),qx(3),0,'Tag','X-axis',vprops{:});
quiver3(o{:},qy(1),qy(2),qy(3),0,'Tag','Y-axis',vprops{:});
quiver3(o{:},qz(1),qz(2),qz(3),0,'Tag','Z-axis',vprops{:});

% Text labels:
tprops = {'Color',clr,'Parent',hl,'HitTest','off','FontSize',12};
offset = scale/10;
tx = origin + qx'*(1+offset/norm(qx));  % The offset just spaces the text
ty = origin + qy'*(1+offset/norm(qy));  %   slightly out past the end of the
tz = origin + qz'*(1+offset/norm(qz));  %   quiver for better readability
text(tx(1),tx(2),tx(3),'x',tprops{:})
text(ty(1),ty(2),ty(3),'y',tprops{:})
text(tz(1),tz(2),tz(3),'z',tprops{:})

end %drawCS()


% ------------------------------------------------------------------------
function [c,p] = roi2cloud(rois)
% Convert roi objects to cloud objects
phases = [rois.Phase];
unique_p = unique(phases);  % also sorts in ascending order
np = numel(unique_p);
for k = np : -1 : 1
    pk = unique_p(k);
    rp = rois(phases == pk);
    % Sort according to slice:
    [~,ind] = sort([rp.Slice]);
    rp = rp(ind);
    % Convert data to 3d & preserve separate traces:
    X = rp(ind).to3dcell;
    % Now create cloud:
    c(k) = Cloud(X);
end
p = unique_p;
end %roi2cloud()


% ------------------------------------------------------------------------
function xyz = rotatePts(xyz,alpha,azel,origin)
%ROTATE Rotate points about specified origin and direction.
%   ROTATE(XYZ,ALPHA,[THETA PHI]) rotates the N-by-3 point cloud XYZ
%   through angle ALPHA about an axis described by the 2-element
%   direction vector [THETA PHI] (spherical coordinates).  
%   All the angles are in radians.  
%
%   THETA is the angle in the xy plane counterclockwise from the
%   positive x axis.  PHI is the elevation of the direction vector
%   from the xy plane (see also SPH2CART).  Positive ALPHA is defined
%   as the righthand-rule angle about the direction vector as it
%   extends from the origin.
%
%   ROTATE(XYZ,ALPHA,[X Y Z]) rotates the points about the direction
%   vector [X Y Z] (Cartesian coordinates). The direction vector
%   is the vector from the center of the plot box to (X,Y,Z).
%
%   ROTATE(...,ORIGIN) uses the point ORIGIN = [x0,y0,y0] as
%   the center of rotation instead the object centroid.
%
%   **Note:
%       This function is partly taken from MATLAB's graphics
%       function ROTATE, but with some modifications for performance &
%       simplicity, and re-ordering of the input arguments (reflected in
%       this help text)

assert(size(xyz,2) == 3, 'XYZ must be n-by-3')

if nargin < 4
    SHIFT = false;
else 
    SHIFT = true;
end

% Find unit vector for axis of rotation
if numel(azel) == 2 % theta, phi
    theta = azel(1);
    phi = azel(2);
    u = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
elseif numel(azel) == 3 % direction vector
    u = azel(:)/norm(azel);
end

% Build rotation matrix:
cosa = cos(alpha);
sina = sin(alpha);
vera = 1 - cosa;
x = u(1);
y = u(2);
z = u(3);
rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';

% Translate if required:
if SHIFT
    xyz = shift(xyz,-origin);
end

% Rotate:
xyz = xyz*rot;

% Translate back, if required:
if SHIFT
    xyz = shift(xyz,origin);
end

end %rotatePts()


% ------------------------------------------------------------------------
function pts = shift(pts,v)
pts = [pts(:,1)+v(1), pts(:,2)+v(2), pts(:,3)+v(3)];
end %shift()


% ------------------------------------------------------------------------
function c = structs2clouds(c,s)
%STRUCTS2CLOUDS Convert structures to clouds (compatability/conversion)
% 
% This function upgrages structs to the current version of the Cloud class.

% Reserve memory:
c = repmat(c,size(s));

% Now we might have spare/old properties that we won't use.  We'll just
% ignore them.
% We'll assign only public or protected properties. Private properties
% should be set some other way, becasue there's a reason that they're
% private (like 'Version', which remains fixed)
mc = metaclass(c);
if isprop(mc,'PropertyList')
    % 2011+ version, a convenient list is provided:
    names = {mc.PropertyList.Name};
    sa = {mc.PropertyList.SetAccess};
else
    % 2010 and earlier version, we must get the stuff manually
    names = getmcprop(mc.Properties,'Name');
    sa = getmcprop(mc.Properties,'SetAccess');
end
    %------------- Legacy helper function (for pre-2011)
    function list = getmcprop(props,field)
        list = cell(1,numel(props));
        for p = 1:numel(props)
            list{p} = props{p}.(field);
        end
    end %getmcprop()

pub  = strcmp(sa,'public');
prot = strcmp(sa,'protected');
names = names( pub | prot );

% Now just process all properties:
for j = 1:numel(names)                  % Cycle through all properties
    if isfield(s,names{j})              % If structure has a field to pull
        [c.(names{j})] = s.(names{j});  % Then pull it
    end
end
end %struct2cloud()

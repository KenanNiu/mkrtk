classdef roi
%ROI class definition.
%
% Constructor Methods:
%  r = roi()                        % Create new 1-by-1 unpopulated roi object 
%  r = roi([])                      % Create an empty 0-by-0 roi object
%  r = roi(s)                       % Create new roi from struct
%  r = roi(s,param,value,...)       % Used for upgrading older structs
%  r = roi(x,y)
%  r = roi(x,y,slice,phase)
%  r = roi(x,y,slice,phase,ImgParent)
%  r = roi(...,PixelSpacing,ImagePositionPatient,ImageOrientationPatient)
%  r = roi(...,ColorSpec)
%  r = roi(...,ColorSpec,LineSpec)
%
% Ordinary Mehtods:
%  r.close          % Close ROI by joining first and last point (if necessary)
%  r.isclosed       % Test to see if roi is closed
%  r.open           % Open a closed ROI (if not alreay open)
%  r.setTag         % Automatically generate a tag
%  
%
% Static Methods:
%  roi.loadfromfile % Load method which will upgrade legacy structs to roi objects
%  roi.loadobj      % Standard load method, called when loading roi objects using LOAD
%  roi.ver          % Return current class version
%
% Getter Methods:
%  r.Area           % Get enclosed area
%  r.Length         % Get length of ROI perimeter


% History
% - Version 2 is the first version to be defined with the class
%   - Version 0 and 1 were structs
% 

% Joshua Martin, May-2012

properties 
    x
    y
    Slice
    Phase
    Image
    Color = 'r';
    LineStyle = '-';
    PixelSpacing
    ImagePositionPatient
    ImageOrientationPatient 
end %properties

properties (Dependent = true, SetAccess = private)
    Area
    Length
end

properties (SetAccess = private)
    Version = 2;        
end %properties (SetAccess = protected)

properties (Transient, SetAccess = private)
    Tag
end %properties (Transient)

methods
    
    %=====================================================================
    % Constructor method:
    
    function r = roi( varargin )    % (constructor)
        if nargin == 0
            % roi()
            return
        elseif nargin == 1 && isempty(varargin{1})
            % roi([]) or roi({})
            r(:) = [];
            return
        elseif nargin >= 1 && isa(varargin{1},'struct')
            % struct
            r = struct2roi( varargin{:});
        elseif nargin >= 2
            % roi(x,y,...)
            [r.x,r.y] = checkxy(varargin{1:2});
            if nargin == 3 || nargin == 6 || nargin == 7 || nargin > 11
                error(message('ROI:roi:IncorectNargin'))
            end
            % roi(x,y,slice,phase)
            if nargin >= 4
                [r.Slice,r.Phase] = checksp(varargin{3:4});
            end
            % roi(x,y,s,p,Image)
            if nargin >= 5
                r.Image = varargin{5};
            end
            % roi(x,y,s,p,img,PixelSpacing,ImagePositionPatient,ImageOrientationPatient)
            if nargin >= 8 
                [r.PixelSpacing, r.ImagePositionPatient,...
                    r.ImageOrientationPatient] = checkpatientcs( varargin{6:8} );
            end
            if nargin >= 9
                r.Color = varargin{9};
            end
            if nargin >= 10
                r.LineStyle = varargin{10};
            end
            r.Tag = genTag();
        else
            error(message('ROI:roi:IncorectNargin'))
        end
    end %roi (constructor)
    
    
    %=====================================================================
    % Ordinary Methods
    
    %----------------------------------------------
    function rset = addpatientcs(rset,dinfo)
        % Insert patient coordinate system data into rois
        ok = isa(dinfo,'struct') &...
            isfield(dinfo,'PixelSpacing') & ...
            isfield(dinfo,'ImagePositionPatient') & ...
            isfield(dinfo,'ImageOrientationPatient');
        if ~ok
            error(message('ROI:addpatientcs:BadDicomInfo'))
        end
        for j = 1:numel(rset)
            dj = dinfo(rset(j).Slice,rset(j).Phase);
            rset(j).PixelSpacing = dj.PixelSpacing;
            rset(j).ImagePositionPatient = dj.ImagePositionPatient;
            rset(j).ImageOrientationPatient = dj.ImageOrientationPatient;
        end        
    end %addpatientcs()
    
    %----------------------------------------------
    function r = close(r)
        if ~isclosed(r)
            TOL = 1e-6;
            gx = diff(r.x([1,end])).^2;
            gy = diff(r.y([1,end])).^2;
            if sqrt(gx.^2 + gy.^2) < TOL    % Almost closed, within tol
                r.x(end) = r.x(1);          % Shift last point to mate
                r.y(end) = r.y(1);          %  with first point
            else                        % Properly not closed
                r.x(end+1) = r.x(1);    % Add endpoint which is 
                r.y(end+1) = r.y(1);    %  same as start point
            end
        end
    end %close()
    
    %----------------------------------------------
    function r = equispace(r,ds)
        % Respample at equal spacing
        if nargin == 1
            ds = 2;
        end
        CL = isclosed(r);
        r = r.open;
        
        dsv = hypot(diff(r.x),diff(r.y));
        s = [0; cumsum(dsv)];
        len = sum(dsv);
        n = ceil( len/ds ) + 1;
        si = linspace(0,len,n)';
        r.x = interp1(s,r.x,si,'pchip'); 
        r.y = interp1(s,r.y,si,'pchip');
        
        if CL
            r = r.close;
        end
        
    end %equispace()
    
    %----------------------------------------------
    function tf = haspatientcs(rset)
        % Check to see if the patient coordinate system has been defined
        ps  = ~cellfun(@isempty,{rset.PixelSpacing});
        ipp = ~cellfun(@isempty,{rset.ImagePositionPatient});
        iop = ~cellfun(@isempty,{rset.ImageOrientationPatient});
        tf  = ps & ipp & iop;
    end %haspatientcs()
    
    %----------------------------------------------
    function r = insert(r,xy)
        % Alias for join
        r = join(r,xy);
    end %insert()
        
    %----------------------------------------------
    function tf = isclosed(r)
        if numel(r.x) <= 2
            tf = false;
        else
            tf = (r.x(1) == r.x(end)) && (r.y(1) == r.y(end));
        end
    end %isclosed()
    
    %----------------------------------------------
    function tf = isclockwise(r)
        tf = curveIsClockwise(r.x,r.y);
    end %isclockwise()
    
    %----------------------------------------------
    function r = join(r,xy)
        % Join one curve segment into roi
            nmin = 3;
        if ~any(size(xy) == 2)
            error(message('ROI:insert:WrongDimensions'))
        end
        
        % Format xy to N-by-2:
        xy = permute(xy, 1+ double(size(xy)==2));
        
        if size(xy,1) < nmin
            error(message('ROI:insert:NotEnoughPoints'))
        end
        
        CL = isclosed(r);
        %r = r.open;
        
        [r.x,r.y] = curveJoin(r.x,r.y,xy(:,1),xy(:,2));
        
        if CL 
            r = close(r);
        end
    end %join()
    
    %----------------------------------------------
    function r = open(r)
        if isclosed(r)
            r.x = r.x(1:end-1);
            r.y = r.y(1:end-1);
        end
    end %open()
    
    %----------------------------------------------
    function r = saveobj(r)
        % If we need to do any special saving manipulation, it can be done
        % here
        
        %if nargin == 1
        %    error(message('ROI:save:FilenameNotSpecified'))
        %end
        %save(filename,'r') %standard matlab save
    end %save()
    
    %----------------------------------------------
    function r = settag(r)
        if isempty(r.Tag)
            r.Tag = genTag();
        end
    end %settag()
    
    %----------------------------------------------
    function xyz = to3d(rset)
        % Convert (x,y) roi(s) to 3D curve using transformation
        %   Input can be an array of roi objects which get merged 
        %   into the one point set for output.
        X = rset.to3dcell;
        xyz = cat(1, X{:});
    end %to3d()
    
    %----------------------------------------------
    function X = to3dcell(rset)
        % Convert (x,y) roi(s) to N-by-1 cell array of (x,y,z) points
        %   Input is normally an array of N roi objects
        %   Output is N-by-1 cell array of (x,y,z) points defining each
        %   input ROI.
        n = numel(rset);
        for j = n : -1 : 1
            r = rset(j);
            % Check integrity:
            [ps,ipp,iop] = checkpatientcs(...
                r.PixelSpacing,r.ImagePositionPatient,r.ImageOrientationPatient);
            % CHECKPATIENTCS will pass if they are all empty, so:
            if isempty(ps) || isempty(ipp) || isempty(iop)
                error(message('ROI:to3d:TransformationNotSpecified'))
            end
            % Check that all transformations in the set are the same:
            if j == n
                PS = ps;
                %IPP = ipp; Image Position Patient differs on each slice
                IOP = iop;
            end
            if ~isequal(ps,PS) || ~isequal(iop,IOP)
                [id,msg] = message('ROI:to3d:DifferentTransformationsInSet');
                warning(id,msg)
            end
            % Convert to 3d and store:
            xy = [r.x(:)'; r.y(:)'] - 1;        %2-by-N, zero-based
            xyzj = xy2xyz(xy, ps, ipp, iop);    %3-by-N
            X{j,1} = xyzj';
        end 
        
    end %to3dcell()
    
    %----------------------------------------------
    function r = toclockwise(r)
        if ~isclockwise(r)
            r.x = flipud(r.x);
            r.y = flipud(r.y);
        end
    end %toclockwise()
    
    %----------------------------------------------
    function r = tocounterclockwise(r)
        if isclockwise(r)
            r.x = flipud(r.x);
            r.y = flipud(r.y);
        end
    end %tocounterclockwise()    
    
    
    %=====================================================================
    % Getter Methods
    
    %----------------------------------------------
    function a = get.Area(r)
        % Calculate area
        xp = r.x;
        yp = r.y;
        a.pixels = polyarea(xp,yp);
        if isempty(r.PixelSpacing)
            a.mm = [];
        else
            xm = xp*r.PixelSpacing(1);
            ym = yp*r.PixelSpacing(2);
            a.mm = polyarea(xm,ym);            
        end
    end %get.Area()
    
    %----------------------------------------------
    function l = get.Length(r)
        % Calculate perimeter length
        dxp = diff(r.x);
        dyp = diff(r.y);
        l.pixels = sum( sqrt(dxp.*dxp + dyp.*dyp) );
        if isempty(r.PixelSpacing)
            l.mm = [];
        else
            dxm = dxp*r.PixelSpacing(1);
            dym = dyp*r.PixelSpacing(2);
            l.mm = sum( sqrt(dxm.*dxm + dym.*dym) );
        end
    end %get.Length()
    
    
end %methods


%=========================================================================
% Static Methods

methods (Static)  
    
    %----------------------------------------------
    function [r,msg] = loadfromfile(fname)
        % Load method wrapper for handling loading of legacy data
        % (pre- version 2.0, before using roi class)
        %
        % MSG can return multiple messages in a cell array indicating if
        % there were problems or cautions with the load/upgrade process.
        if exist(fname,'file') ~= 2
            error(message('ROI:loadfromfile:CouldNotReadFile',fname))
        end
        msg = {};
        msgData = {};
        data = load(fname);
        % Now we concatenate all roi objects onto the output, and
        % upgrade all structs with a warning (if necessary)
        r = roi([]);
        fields = fieldnames(data);
        for j = 1:numel(fields)
            fieldj = data.(fields{j});
            if isa(fieldj,'roi')
                % Got a roi which is fully qualified
                n = numel(fieldj);
                r(end+1:end+n) = data.(fields{j});
            elseif isa(fieldj,'struct') && isfield(fieldj,'Slice')
                % Got a structure to convert/upgrade
                if isfield(fieldj,'Phase')
                    tmp = roi(fieldj);
                    n = numel(tmp);
                    r(end+1:end+n) = tmp;
                else
                    msg{end+1} = ['Some ROIs have not been loaded because ',...
                        'they do not have the "Phase" paremeter specified. ',...
                        'Upgrade the files and try again.']; %#ok<AGROW>
                    msgData{end+1} = 'Phase'; %#ok<AGROW>
                end
            end
        end %for
        % Now add warning messages for problems that have been encountered:
        if any(cellfun(@isempty,{r.PixelSpacing})) || ...
                any(cellfun(@isempty,{r.ImagePositionPatient})) || ...
                any(cellfun(@isempty,{r.ImageOrientationPatient}))
            msg{end+1} = sprintf(['Some ROIs have been loaded and do not have',...
                ' one or more of the following properties specified:',...
                '\n  PixelSpacing',...
                '\n  ImagePositionPatient',...
                '\n  ImageOrientationPatient',...
                '\nThese will need to be specified before converting to 3D.']);
            msgData(end+1:end+3) = {'PixelSpacing' 'ImagePositionPatient' 'ImageOrientationPatient'};
        end
        % Now if we have a problem, and the request has been made without
        % calling for the error message, we throw an error, because it
        % would otherwise go unnoticed: 
        if ~isempty(msg) && nargout < 2
            error(message('ROI:struct2roi:DataRequired',msgData{:}))
        end
    end %loadfromfile()
    
    %----------------------------------------------
    function r = loadobj(r)
        % This method is called by default when any object of this class is
        % loaded.  If the automatic load succeeds, we have successfully
        % created the object.  If it fails, we will have a struct
        if isa(r,'struct')
            % Automatic load failed.
            keyboard
            % Need a try/ catch to inform the user that they might need to
            % upgrade with parameters?
            r = struct2roi(r);
        end
        r.Tag = genTag();
        %r = loadFromFile(fname);
    end %loadobj()
    
    %----------------------------------------------
    function pt3d = pix2mm3D(pt,PixelSpacing,ImagePositionPatient,ImageOrientationPatient)
        % This static method gives direct access to the xy2xyz function
        % which converts 2d image points to 3d measurements
        pt3d = xy2xyz(pt,PixelSpacing,ImagePositionPatient,ImageOrientationPatient);
    end %pix2mm3D()
    
    %----------------------------------------------
    function v = ver()
        % This static method gives direct access to the version property to
        % enable checking of the current class version
        v = roi.Version;
    end %ver()
    
end %methods(Static)


end %classdef roi


% ========================================================================
% Helper Functions


% ------------------------------------------------------------------------
function [ps, ipp, iop] = checkpatientcs(ps,ipp,iop)

if isempty(ps) && isempty(ipp) && isempty(iop)
    % ok
    return
elseif isempty(ps) || isempty(ipp) || isempty(iop) 
    error(message('ROI:checkpatientcs:RequirePsIppIop'))
end

if numel(ps) ~= 2 || ~isa(ps,'float')
    error(message('ROI:checkpatientcs:IncorrectPixelSpacing'))
end

if numel(ipp) ~= 3 || ~isa(ipp,'float')
    error(message('ROI:checkpatientcs:IncorrectImagePositionPatient'));
end

if numel(iop) ~= 6 || ~isa(iop,'float')
    error(message('ROI:checkpatientcs:IncorrectImageOrientationPatient'))
end

% More easily readable as 1-by-N:
ps  = ps(:);
ipp = ipp(:);
iop = iop(:);

end %checkpatientcs()


% ------------------------------------------------------------------------
function [slice,phase] = checksp(slice,phase)
% Check integrity of slice and phase numbers
if ( isempty(slice) && isempty(phase) ) 
    %ok
    return
end
if mod(slice,1) ~= 0
    [id,msg] = message('ROI:checksp:NonIntegerValue');
    warning(id,msg)
    slice = round(slice);
end
if mod(phase,1) ~= 0
    [id,msg] = message('ROI:checksp:NonIntegerValue');
    warning(id,msg)
    phase = round(phase);
end

end %checksp()
    

% ------------------------------------------------------------------------
function [x,y] = checkxy(x,y)
% Ensure x & y are the same size:
if numel(x) ~= numel(y)
    error(message('ROI:checkxy:InconsistentLengths'))
end

% Ensure x is a vector:
if length(find(size(x)>1))>1
    error(message('ROI:checkxy:XNotVector'))
end

% Ensure minimum number of points to define an roi:
if numel(x) < 3
    error(message('ROI:checkxy:NotEnoughPoints'))
end

% Format output
x = x(:);
y = y(:);

end %checkxy()


% ------------------------------------------------------------------------
function tf = curveIsClockwise(x,y)
mx = mean(x);
my = mean(y);
tf = sum( unwrap( diff( angle( complex(x-mx, y-my) ) ) ) ) < 0;
end %curveIsClockwise()


% ------------------------------------------------------------------------
function [newx,newy] = curveJoin(x,y,xs,ys)
% Joint segment (xs,ys) into the curve (x,y)

% Enforce column vectors:
x = x(:);
y = y(:);
xs = xs(:);
ys = ys(:);

% Format as clockwise:
CW = curveIsClockwise(x,y);
if ~CW
    x = flipud(x);
    y = flipud(y);
end
if ~curveIsClockwise(xs,ys)
    xs = flipud(xs);
    ys = flipud(ys);
end

% Record first point:
xfirst = x(1);
yfirst = y(1);

% Find intersections:
[~,~,ia,ib] = intersectCurves(x,y,xs,ys);

% Now trim / extend the segment to the intersection points:
np = round(abs(diff(ib)));               % Required number of points:
t = linspace(ib(1),ib(end),np);     
xs = interp1(xs,t,'linear','extrap');
ys = interp1(ys,t,'linear','extrap');

% Two segments of the main curve:
seg1x = x(ceil(ia(1)):floor(ia(end)));
seg1y = y(ceil(ia(1)):floor(ia(end)));

seg2x = [x(ceil(ia(end)):end); x(1:floor(ia(1)))];
seg2y = [y(ceil(ia(end)):end); y(1:floor(ia(1)))];

% Join the insertion segment onto the segments of the main curve:
[seg1x,seg1y] = joinsegs(seg1x,seg1y,xs(:),ys(:));
[seg2x,seg2y] = joinsegs(seg2x,seg2y,xs(:),ys(:));

% Result should be the largest of the two areas:
if polyarea(seg1x,seg1y) > polyarea(seg2x,seg2y)
   newx = seg1x;
   newy = seg1y;
else
   newx = seg2x;
   newy = seg2y;
end

% Restore first point:
%   Curve has been resampled, so it won't match up exactly.
%   Do a min distance search:
[~,i0] = min( hypot(newx-xfirst,newy-yfirst) );
resort = [i0:numel(newx), 1:i0-1];
newx = newx(resort);
newy = newy(resort);

% Revert to CCW if that was what was input
if ~CW
    newx = flipud(newx);
    newy = flipud(newy);
end


    %------------------------------
    function [xnew,ynew] = joinsegs(xmain,ymain,xseg,yseg)
        d1 = hypot(xmain(1)-xseg(end),ymain(1)-yseg(end) );
        d2 = hypot(xmain(1)-xseg(1),ymain(1)-yseg(1));
        if d1 < d2  % curve(1) near segment(end)
            xnew = [xmain; xseg(:)];
            ynew = [ymain; yseg(:)];
        else        % curve(end) near segment(1)
            xnew = [xmain; flipud(xseg(:))];
            ynew = [ymain; flipud(yseg(:))];
        end
        
        % Remove duplicate points (if there are any)
        xy = [xnew,ynew];
        %xy = unique(xy,'rows','stable');       % This only came in in 2012
        [~,i1] = unique(xy,'rows');             % Instead we do this...
        xy = xy(sort(i1),:);                    %  ...then restore the original order...
                
        xnew = xy(:,1);
        ynew = xy(:,2);
                
    end %joinsegs()

end %curveJoin()


% ------------------------------------------------------------------------
function tag = genTag()
persistent tagnum
if isempty(tagnum)
    tagnum = 1;
end
tag = sprintf('ROI-%d',tagnum);
tagnum = tagnum+1;
end %genTag()


% ------------------------------------------------------------------------
function [xi,yi,ia,ib] = intersectCurves(xa,ya,xb,yb)
% Find curve segment intersections.
%   ia - fractional index of intersection position along curve (xa,ya)
%   ib - fractional index of intersection position along curve (xb,yb)
%
% More extensive & robust procedures can be found in polyxpoly or on the
% file exchange in submission #11837, but we have some differences here

xa = xa(:);
ya = ya(:);
xb = xb(:);
yb = yb(:);

[xi(1),yi(1),ia(1),ib(1)] = getIntersection(xa,ya,xb,yb,'first');
[xi(2),yi(2),ia(2),ib(2)] = getIntersection(xa,ya,xb,yb,'last');

% Sort in terms of their position along curve a:
[ia,ind] = sort(ia);
xi = xi(ind);
yi = yi(ind);
ib = ib(ind);

    % -------------------------------------------------------------
    function [xi,yi,ia,ib] = getIntersection(x0,y0,x1,y1,type)
        
        if isequal(type,'last')
            x1 = flipud(x1);
            y1 = flipud(y1);
            n1 = numel(x1);
        end
        % Only look at the first part of the curve:
        SPAN = 0.2; % percentage of the curve to examine 
        ne = ceil(SPAN*numel(x1))+1;
        x1 = x1(1:ne);
        y1 = y1(1:ne);
        
        % We're looking for the first solution we can get, in the following
        %   order of priority:
        % First priority:
        [xi,yi,ia,ib] = getExplicitIntersection(x0,y0,x1,y1,true);
        % Second priority:
        %if isempty(xi)
        %    [xi,yi,ia,ib] = getExplicitIntersection(x0,y0,x1(1:2),y1(1:2),false);
        %end
        % Third priority:
        if isempty(xi)
            [xi,yi,ia,ib] = getNearSolution(x0,y0,x1,y1);
        end
        
        % Cleanup
        if isequal(type,'last')
            ib = n1+1-ib;
        end
        
    end %getIntersection()
        
    % -------------------------------------------------------------
    function [xi,yi,ia,ib] = getExplicitIntersection(xa,ya,xb,yb,explicit)
        % Get explicit crossings of the two curves
        [xi,yi,ia,ib] = deal([]);
        for j = 1:numel(xa)-1
            for k = 1:numel(xb)-1
                x1 = [xa(j), xa(j+1)];
                y1 = [ya(j), ya(j+1)];
                x2 = [xb(k), xb(k+1)];
                y2 = [yb(k), yb(k+1)];
                
                xk = [ x1(:) x2(:) ];
                yk = [ y1(:) y2(:) ];
                
                dx = diff(xk,1);
                dy = diff(yk,1);
                
                den = dx(1)*dy(2)-dy(1)*dx(2);
                if abs(den) < eps
                    continue
                end
                ua = (dx(2)*(yk(1)-yk(3))-dy(2)*(xk(1)-xk(3)))/den;
                ub = (dx(1)*(yk(1)-yk(3))-dy(1)*(xk(1)-xk(3)))/den;
                % Check if intersection lies within segments:
                isInSegment = all(([ua ub] >= 0) & ([ua ub] <= 1));
                if (explicit && isInSegment) || ~explicit
                    % ~explicit will be used in the case where line a
                    % contains only one segment
                    xi = xk(1)+ua*dx(1); 
                    yi = yk(1)+ua*dy(1);
                    ia = j + ua;
                    ib = k + ub;
                    % The first solution is all we want; bail when we've
                    % got it:
                    return
                end
                
            end
        end
    end %getExplicitIntersection()


    % -------------------------------------------------------------
    function [xi,yi,ia,ib] = getNearSolution(x0,y0,x1,y1)
        dmin = inf;
        for k = 1:numel(x1)
            dv = hypot(x1(k)-x0,y1(k)-y0);
            d = min(dv);
            if d < dmin
                j = find(dv==d,1,'first');
                xi = x0(j);
                yi = y0(j);
                ia = j;
                ib = k;
                dmin = d;
            end
        end            
    end %getNearSolution()


end %intersectCurves()


% ------------------------------------------------------------------------
function s = mergeStructs(a,b)
% Values in B overwrite values in A, unless empty
bfields = fieldnames(b);

% Inititalise:
s = a;

% Insert fields if they don't exist, or if b's field is non-empty:
for j = 1:numel(bfields)
    if ~isfield(a,bfields{j}) || ~isempty(b.(bfields{j}))
        s.(bfields{j}) = b.(bfields{j});
    end
end
end %mergeStructs()


% ------------------------------------------------------------------------
function varargout = message(msgid,varargin)
% Generate warning or error messaged from message id
switch msgid
            
    case 'ROI:addpatientcs:BadDicomInfo'
        message = ['Second input should be a struct containing the fields ',...
            'PixelSpacing, ImagePositionPatient, and ImageOrientationPatient.'];
        
    case 'ROI:checkxy:NotEnoughPoints'
        message = 'At least 3 points are required to define a ROI';
        
    case 'ROI:checkpatientcs:IncorrectImagePositionPatient'
        message = 'ImagePositionPatient must be a 3-by-1 float.';
        
    case 'ROI:checkpatientcs:IncorrectImageOrientationPatient'
        message = 'ImageOrientationPatient must be a 6-by-1 float.';
        
    case 'ROI:checkpatientcs:IncorrectPixelSpacing'
        message = 'PixelSpacing must be a 2-element float.';
        
    case 'ROI:checkpatientcs:RequirePsIppIop'
        message = 'PixelSpacing, ImagePositionPatient and ImageOrientationPatient must be provided together';
        
    case 'ROI:checksp:NonIntegerValue'
        message = 'Slice and phase indices should be integers. Value rounded';
    
    case 'ROI:checkxy:InconsistentLengths'
        message = 'X and Y must be the same length.';
        
    case 'ROI:checkxy:XNotVector'
        message = 'X and Y must be vectors';
        
    case 'ROI:insert:NotEnoughPoints'
        message = 'At least 3 points are required to modify a ROI.';
        
	case 'ROI:insert:WrongDimensions'
        message = 'XY must be N-by-2 or 2-by-N set of points.';
        
    %case 'ROI:load:RequireMatFile'
    %    message = 'Filename of a MAT file must be provided.';
        
    case 'ROI:loadfromfile:CouldNotReadFile'
        message = sprintf('Unable to read file %s: No such file or directory.',varargin{1});
        
    case 'ROI:MergeStructs:BadOpt'
        message = sprintf('Bad merge option: %s.',varargin{1});
        
    case 'ROI:roi:IncorectNargin'
        message = 'Incorrect number of input arguments.';
        
    case 'ROI:roi:NdStructArraysNotSupported'
        message = ['Creation of ROIs from N-dimensional structures are not yet supported.',...
            ' Create each ROI with an individual call to the class (constructor).'];
        
    case 'ROI:save:FilenameNotSpecified'
        message = 'Filename must be specified.';
        
    case 'ROI:struct2roi:DataRequired'
        message = ['You must specify the parameter(s) ',sprintf('"%s", ',varargin{:}),...
            'in order to upgrade this ROI definition to the current version.'];
        
    case 'ROI:struct2roi:UnhandledVersion'
        message = sprintf('Unhandled version [%d].',varargin{1});
        
    case 'ROI:to3d:SingleObjectOnly'
        message = 'The "to3d" method works only on one object at a time.';
        
    case 'ROI:to3d:DifferentTransformationsInSet'
        message = ['The ROIs in the set specified for conversion to 3D ',...
            'had different transformations. Normally this is undesirable ',...
            'since the result is merged into a single point set.'];
        
    case 'ROI:to3d:TransformationNotSpecified'
        message = ['Transformation parmeters not specified.  Please specify ',...
            'PixelSpacing, ImagePositionPatient, and ImageOrientationPatient',...
            ' before trying to convert to 3D.'];
        
    otherwise
        message = sprintf('Unhandled Message identifier: %s',msgid);
end

switch nargout
    case 1
        % Error message:
        msg.message = message;
        msg.identifier = msgid;
        varargout = {msg};
        
    case 2
        % Warning message:
        varargout = {msgid,message};
end

end %message()


% ------------------------------------------------------------------------
function s = ratifyStruct(s,VSN,varargin)
% Additional parameters may be required to qualify a struct to the current
% version (VSN).  If this is the case 
% required and not provided, an error will be thrown.  The additional data
% is provided with parameter / value pairs, as in the example that follows:


% These are the only inputs which we will handles.  Others may be provided,
% but they will be ignored by the parser:
p = inputParser;
addParamValue(p,'Phase',[],@isnumeric)
addParamValue(p,'Image',[],@ischar)
addParamValue(p,'ImagePositionPatient',[],@(x) isnumeric(x));
addParamValue(p,'ImageOrientationPatient',[],@(x) isnumeric(x))
addParamValue(p,'PixelSpacing',[],@(x) isnumeric(x));
parse(p,varargin{:})    % Run the parser
Res = p.Results;        % Take a copy so we can do our own parsing & modify

% Run checks on data:
[Res.PixelSpacing, Res.ImagePositionPatient, Res.ImageOrientationPatient] = ...
    checkpatientcs(Res.PixelSpacing, Res.ImagePositionPatient, Res.ImageOrientationPatient);

% Should these checks (including checkixp & checkps) be moved into setter
% methods?

v = 0;
if isfield(s,'Version')
    v = s.Version;
elseif isa(s,'roi')
    v = s.Version;
end


% Manage all possible combinations of upgrades:
if v == 0 && VSN == 2       % 0 -> 2
    s = ver0to1(s);
    s = ver1to2(s);
elseif v == 1 && VSN == 2   % 1 -> 2
    s = ver1to2(s);
elseif v == 2 && VSN == 2   % 2 -> 2
    s = mergeStructs(sVersionDef(2,[]),s); 
else
    error(message('ROI:struct2roi:UnhandledVersion',VSN))
end

% Merge input data, if provided:
s = mergeStructs(s,Res);

% Now we have a structure which contains all the fields, and perhaps some
% old or redundant ones.

    % -------------------------------------------   
    function s = ver0to1(s)
        % Upgrade from version 0 to version 1
        % Going from 0 -> 1, we added:
        %   DcmFile
        %   Version
        %
        %s.DcmFile = [];
        s = mergeStructs(sVersionDef(1,[]),s);
        s.Version = 1;
    end %ver0to1()

    % -------------------------------------------
    function s = ver1to2(s)
        % Upgrage from version 1 to version 2
        %
        % Unfortunately, there was a change without a version increment, so
        % that if we have a version 1 ROI, it may or may not contain:
        %   Phase
        %
        % If Phase exists, use it.  If not, require it:
        if ~isfield(s,'Phase') || isempty(s.Phase)
            if isempty(Res.Phase)
                error(message('ROI:struct2roi:DataRequired','Phase'))
            else
                s.Phase = Res.Phase;
            end
        end

        % In version 2, 'DcmFile' becomes 'Image':
        s.Image = s.DcmFile;
        s = rmfield(s,'DcmFile');
        
        % 'Interp' was removed:
        %s = rmfield(s,'Interp');
        
        % And we have added the following fields:
        %   ImagePositionPatient
        %   ImageOrientationPatient
        %   PixelSpacing
        s = mergeStructs(sVersionDef(2,[]),s);
        s.Version = 2;
    end %ver1to2()

end %qualifyStruct()


% ------------------------------------------------------------------------
function r = struct2roi(s,varargin)
%STRUCT2ROIS Compatability/conversion function
%
% This function can upgrage structs to the current version of the roi
% class.  In order to upgrade older structs, some additional info may be
% required.  See qualifyStruct for what that info is.
%
% See also Cloud.m>structs2clouds for a similar implementation

% Instigate empty rois:
r = repmat(roi(),size(s));

% Now we might have spare/old properties that we won't use.  We'll just
% ignore them.
% We'll assign only public or protected properties. Private properties
% should be set some other way, becasue there's a reason that they're
% private (like 'Version', which remains fixed)
mc = metaclass(r);
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

% Now process all elements of the structure & turn them into roi objects:
for k = 1:numel(s)
    r(k).Tag = genTag();
    
    % First we need to make the structure a fully ratified according to the
    % current version of the class.  If ratifyStruct is not up to date with the
    % current version, it will throw an error.
    tmp = ratifyStruct(s(k),r(k).Version,varargin{:});
    
    % Then import the properties:
    for j = 1:numel(names)
        r(k).(names{j}) = tmp.(names{j});
    end
    
end

% Done!    
end %struct2roi()


% ------------------------------------------------------------------------
function s = sVersionDef(VSN,fill)
% Structure Version Definition
%
% Here we define what the structure looked like for each version
%
% if fill = {}, Returns an empty 0x0 structure
% if fill = [], Returns a 1x1 structure with all fields [];
if nargin == 1;
    fill = {};
end
switch VSN
    case 0
        fields = {...
            'x';
            'y';
            'Slice';
            'Interp';
            'Color';
            'LineStyle';
            'Tag'};
            
    case 1
        fields = {...
            'x';
            'y';
            'Slice';        
            'Phase';        % ADDED (may or may not exist in a Version 1 roi struct)
            'DcmFile';      % ADDED
            'Interp';
            'Tag';
            'Color';
            'LineStyle';
            'Version'};
        
    case 2
        fields = {...
            'x';
            'y';
            'Slice';
            'Phase';
            'Image';
            'Color';
            'LineStyle';
            'PixelSpacing';             % ADDED
            'ImagePositionPatient';     % ADDED
            'ImageOrientationPatient';  % ADDED
            'Version';
            'Tag'};
        % Dependent Properties (not saved or restored to/from file):
        %   'Area' 
        %   'Length'
        
        % Removed properties:
        %   'Interp'
        
    otherwise
        error('Unandled Version')
end
sargs = fields(:)';
sargs(2,:) = {fill};
s = struct(sargs{:});
end %sVersionDef()


% ------------------------------------------------------------------------
function [xyz] = xy2xyz(pt,PixelSpacing,ImagePositionPatient,ImageOrientationPatient)
%XY2XYZ Convert image pixel locations to 3D locations in the patient
% coordinate system.
% 
% INPUTS:
%                        pt: 2-by-N list of (i,j) pixel locations, ZERO-BASED
%              PixelSpacing: 2-by-1 property from DICOM header
%      ImagePositionPatient: 3-by-1 property from DICOM header
%   ImageOrientationPatient: 6-by-1 property from DICOM header
%
%
% OUTPUT:
%           xyz: 3-by-N list of (x,y,z) mm locations in the patient coordinate system
%
% This function essentially does the homogeneous matrix multiplication of
% the form:
%
%  +- -+          +- -+
%  | x |    +- -+ | i |
%  | y | =  | M | | j |
%  | z |    +- -+ | 0 |
%  | 1 |          | 1 |
%  +- -+          +- -+
%
% The input image point PT is described as (i,j) and according to the
% standard is a zero-based index.  Hence:
%
%   - (0,0) is top left voxel
%   - ImagePositionPatient describes the centre of top-left voxel
%
%
% The details are described in the DICOM standard, Section C.7.6.2.1.1,
% page 275: 
% http://medical.nema.org/dicom/2004/04_03PU.PDF


% Check shape of input points:
assert(size(pt,1)==2, 'Image points must be provided in 2-by-N form of (x,y) pairs')

% Location of centre of top-left voxel:
Sx = ImagePositionPatient(1);
Sy = ImagePositionPatient(2);
Sz = ImagePositionPatient(3);

% Row direction cosine:
Xx = ImageOrientationPatient(1);
Xy = ImageOrientationPatient(2);
Xz = ImageOrientationPatient(3);

% Column direction cosine:
Yx = ImageOrientationPatient(4);
Yy = ImageOrientationPatient(5);
Yz = ImageOrientationPatient(6);

% PixelSpacing:
di = PixelSpacing(1);
dj = PixelSpacing(2);

% Mapping matrix:
M = [...
    Xx*di Yy*dj 0 Sx
    Xy*di Yy*dj 0 Sy
    Xz*di Yz*dj 0 Sz
      0     0   0  1 ];

% Compute the matrix multiplication
%
%   Basic form is as follows:  
%       Pxyz = M*[pt(:);0;1];
%       xyz = Pxyz(1:3);
%
% In vectorised form:
Pij =  zeros(4,size(pt,2));
Pij(1:2,:) = pt;
Pij(4,:)   = 1;

Pxyz = M*Pij;

xyz = Pxyz(1:3,:);
end %xy2xyz()


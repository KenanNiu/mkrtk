function M = curve_midline(X)
% X is an N-by-1 cell array of M-by-3 sets of points which represent
% different curves
%
if iscell(X)
    % Recursion:
    M = cell(size(X));
    for j = 1:numel(X)
        print_status()
        M{j} = curve_midline(X{j});
    end
    
elseif size(X,1) < 4
    M = X;
    
else
    
    % We have to find the midline in 2D
    p = planeThroughTrace(X);
    
    % And the rotation matrix:
    R(:,1) = p(4:6);
    R(:,2) = p(7:9);
    R(:,3) = cross(p(4:6), p(7:9));
    
    % Transform 3D points on a constant z-plane:
    XPlanar = X*R;
    X2 = XPlanar(:,1:2);
    
    % Get midline of 2d points:
    M2 = skel_midline2d(X2);
    %M2 = voronoi_midline2d(X2);
    %M2 = delaunay_midline2d(X2);
            
    % Transform back to 3D:
    M3 = M2;
    M3(:,3) = XPlanar(1,3);
    M3 = M3*R';
    
    % And that's our result:
    M = M3;
    
    %{
    %% Plot
    figure, plot3(X(:,1),X(:,2),X(:,3),'k.','MarkerSize',7), grid on, axis equal
    hold on
    plot3(M3(:,1),M3(:,2),M3(:,3),'g-','LineWidth',2)
    view(90,0)
    %}
    
end
    % ---------------------
    function print_status()
        ds = dbstack;
        nt = sum(strcmp({ds.name},mfilename))-1;
        indent = repmat('  ',1,nt);
        fprintf('%sCalculating midline %d of %d\n',indent,j,numel(X));
    end

end %curve_midline()

% ------------------------------------------------------------------------
function m = skel_midline2d(xy)
% first we create an representative image which describes the polygon:

[B,T,S] = poly2bw(xy);

% Internal distance map:
%D = bwdist(~B);

% Skeletonize:
BW3 = bwmorph(B,'skel',Inf);

% We would probably do well to trim some of the branches off
% Convert back to curve:
np = size(B,2)/20;
[yi,xi] = find(BW3);                            % Convert to curve
yi = smooth(xi,yi,np,'rlowess');      % Smooth
m = [xi yi]*S + repmat(T(1,:),numel(xi),1);     % Re-scale back to original

%{
figure, 
imshow(B)

figure,
imshow(BW3);
hold on
plot(xi,yi,'g-')
%%
figure
plot(xy(:,1),xy(:,2),'b.','MarkerSize',7);
hold on
plot(m(:,1),m(:,2),'r-','Marker','x')
%}

end %()


% ------------------------------------------------------------------------
function [BW,T,S] = poly2bw(xy)
%
%
% Image points can be returned to their original poly space by doing:
%   xy = [r c]*S + T(r,:)
%%
x = xy(:,1);
y = xy(:,2);
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

% Maximum dimension of 1000px:
N = 500;
rx = range(x);
ry = range(y);
if rx > ry
    nx = N;
    ny = round( nx *ry/rx );
else
    ny = N;
    nx = round( ny *rx/ry );
end

% Scale
sx = (xmax-xmin)/nx;
sy = (ymax-ymin)/ny;
S = diag([sx sy]);

% Shift:
T = zeros(size(xy));
T(:,1) = xmin;
T(:,2) = ymin;


% Apply transformations to convert to image space:
xyi = (xy - T)*(inv(S));

% Mask:
BW = poly2mask(xyi(:,1),xyi(:,2),ny,nx);

%{
figure, imshow(BW,[]), hold on, plot(xyi(:,1),xyi(:,2),'g'), set(gca,'Visible','on')
%%

figure, 
ht = plot(xy(:,1),xy(:,2),'r*');
axis auto, axis equal, drawnow, axis manual
pause
for j = 1:size(xy,1)
    set(ht,'XData',xy(1:j,1),'YData',xy(1:j,2))
    pause(0.05)
end
%}
%%
if all(~BW(:))
    fprintf(2,'poly2mask failed to produce a result\n')
    keyboard
end

end

% ------------------------------------------------------------------------
function m = delaunay_midline2d(xy)
x = xy(:,1);
y = xy(:,2);

k = [1:numel(x), 1]';
Constraints = [k(1:end-1) k(2:end)];
dt = DelaunayTri(x,y,Constraints);
inside = dt.inOutStatus();

tr = TriRep(dt(inside, :), dt.X);

figure
triplot(tr,'r');
grid on, axis tight, hold on

% Circumcentres lie on the medial axis:
cc = tr.circumcenters();

% Now split every line in every tri:
mp = [];
cpts = x + 1i*y;
for j = 1:size(tr.Triangulation,1)
    cptsj = cpts(tr.Triangulation(j,:));
    l = [...
        abs(diff(cptsj([1,2])))
        abs(diff(cptsj([2,3])))
        abs(diff(cptsj([3,1])))];
    [~,idx] = sort(l,'descend');
    mpj = [...
        mean(cpts(tr.Triangulation(j,[1,2])))
        mean(cpts(tr.Triangulation(j,[2,3])))
        mean(cpts(tr.Triangulation(j,[3,1])))];
    mp(end+1:end+2,:) = mpj(idx(1:2));    
end

error('This has major problems with traces that have two sections')

% Now do a robust smooth on combined data:
% vxy = complex(vx(:),vy(:));
% C = [mp; vxy];
% C = unique(C,'rows');
% 
% z = smoothn(C,'robust');
% 
% plot(xy(:,1),xy(:,2),'k.','MarkerSize',7), grid on, axis equal, hold on
% plot(real(C),imag(C),'bo')
% 
% plot(real(z),imag(z),'g','LineWidth',2)
% 
% m = [real(z), imag(z)];


end



% ------------------------------------------------------------------------
function m = voronoi_midline2d(xy)
% Voronoi midline uses a voronoi diagram to ascertain the midline of the
% curve defined by the N-by-2 point set xy (normally a tenondon tracing).
%
% If the curve defined by xy has actually two subcurves, invariably only
% the larger of the two will be recogonised by the algorithm as the graph
% trimming algorithms tend to split the two components apart.


%profile on

%

% Voronoi diagram has the midline in there somewhere:
[vx,vy] = voronoi(xy(:,1),xy(:,2));

% Remove NaNs, Infs, and any other things that voronoi occasionally
% creates:
ok = all( isfinite(vx) & isfinite(vy), 1);
vx = vx(:,ok);
vy = vy(:,ok);

% And the really long lines aren't what we want, so chuck them too
l = sqrt( diff(vx).^2 + diff(vy).^2 );
ok = l < (mean(l) + 3*std(l));
vx = vx(:,ok);
vy = vy(:,ok);

%{
%% Get ourselves a convex hull so that we have a closed boundary so we
% can also throw away any voronoi points outside the hull:
k = convhull(xy);
K = xy(k,:);
in = inpolygon(vx,vy,K(:,1),K(:,2));
in = all(in,1);
vx = vx(:,in);
vy = vy(:,in);
%}

%
% Throw away points outside the curve
in = inpolygon(vx,vy,xy(:,1),xy(:,2));
in = all(in,1);
vx = vx(:,in);
vy = vy(:,in);
%}

% Trim branches, leaving the trunk
[vx,vy] = erode_voronoi_branches(vx,vy);          % This is fast, but still leaves problems
z = explode_voronoi_branches(complex(vx,vy));     % But this is slow to use on its own

% Check ends:
%z = trim_cooked_ends(z);

%toc
%profile off
%profile viewer

% Convert to a line
id1 = points_with_occurrence(real(z),imag(z),@(n)n==1);
z = trace_from_end(z,id1(1));

% Now do a robust smooth:
z = smoothn(z,'robust');

% This will do for now:
m = [real(z), imag(z)];

%{
% 2D plot:
figure(42), plot(xy(:,1),xy(:,2),'k.','MarkerSize',7), grid on, axis equal, hold on
plot(vx(1,:),vy(1,:),'r*')
lims = axis(gca);
plot(vx,vy)
axis(gca,lims)
plot(real(z),imag(z),'g-','LineWidth',2)
hold off
keyboard
%}


end %midline2d()


% % ------------------------------------------------------------------------
% function z2 = segments2path(z)
% % Convert un-ordered voronoi segments to a full path.
% 
% % Check it's ready:
% id1 = points_with_occurrence(real(z),imag(z),@(n)n==1);
% idn = points_with_occurrence(real(z),imag(z),@(n)n>2);
% assert(numel(id1) == 2, 'Voronoi segment set has more than two branch ends')
% assert(numel(idn) == 0, 'Voronoi segment set has one or more unhandled junctions')
% 
% % Trace out:
% n = size(z,2)+1;
% Fr
% 
% 
% end

% ------------------------------------------------------------------------
function z = explode_voronoi_branches(z)
% This explodes off all branches at 3-connected sites or higher
% Is recursive and very slow...

% As a first step, remove any non-connected branches by selecting only the
% largest component:
lbl = vdlabel(z);

% Trim the voronoi diagram to keep only the dominant component:
keepid = mode(lbl);
z(:,lbl ~= keepid) = [];
lbl(lbl ~= keepid) = [];

% Safety check:
if ~isempty( points_with_occurrence(real(z),imag(z),@(n)n>=4) )
    fprintf(2,'This is a special case which should be handled, but we must check...\n');
    keyboard
end

% Find all forks:
fids = points_with_occurrence(real(z),imag(z),@(n)n>=3);
% Locations of forks:
fz = unique(z(fids));

for j = 1:numel(fz)
    % At each fork:
    %   - progress through each joined segment, and 
    %       - remove them one at a time
    %       - do a vdlabel and label all connected components
    %       - do a length summation of each labelled segment
    %       - store z definition which includes only the largest component
    %   - Then compare lengths of all produced voronoi diagrams & choose largest
    
    % Find columns of all segments joined at this junction:
    [~,c] = find(z==fz(j));
    if isempty(c)
        % Components have already been removed by a previous trim
        continue
    end
    zj = cell(1,numel(c));
    trunk_size = zeros(1,numel(c));
    for k = 1:numel(c)
        % Drop each branch in turn:
        zk = z;
        zk(:,c(k)) = [];
        
        % Label the result:
        lblk = vdlabel(zk);
        % Select just the trunk:
        zk = zk(:,lblk==mode(lblk));
        
        % Store:
        zj{k} = zk;
        trunk_size(k) = vdlength(zk);
        
    end
    % Select the longest of the resulting diagrams:
    [~,id] = max(trunk_size);
    
    %{
    figure(99)
    plot(z,'Color',[1 1 1]*0.7,'LineWidth',1.8)
    hold on
    plot(zj{id},'-','Color','r')
    title('vdlabel was producing bad results and dropping parts')
    %}
    
    %keyboard
    z = zj{id};
    
    %figure(99)
    %plot(z,'Color','r')
    %pause(0.1)
end

end


% ------------------------------------------------------------------------
function s = vdlength(z)
% length of all segments in complex voronoi diagram
s = sum(abs(diff(z)));

end


% ------------------------------------------------------------------------
function labels = vdlabel(z)
% Single pass component labelling, adapted from:
%   http://en.wikipedia.org/wiki/Connected-component_labeling#One-pass_version

%zin = z;
N=size(z,2);
Connected = zeros(1,N);
Mark = 1;
Difference = 1;
No_of_Objects = 0;
sz = size(z);
F = false(size(z));

for j = 1:N
    No_of_Objects = No_of_Objects +1;
    Index = j;
    if Connected(Index) == 0
        Connected(Index) = Mark;
    end
    
    while ~isempty(Index)
        % Current point:
        pts = z(:,Index);
        z(:,Index) = NaN;
        
        %{
        % Neighbours:
        neigh_cols = [];
        for k = 1:numel(Index)
            [~,nck] = find( z==pts(1,k) | z==pts(2,k) );
            neigh_cols(end+1:end+numel(nck)) = nck;
        end
        Neighbours = unique(neigh_cols);
        Index = Neighbours;
        Connected(Index)=Mark;
        %}
        
        % Neighbours:
        neigh_tf = F;
        for k = 1:numel(Index)
            neigh_tf = ( neigh_tf | z==pts(1,k) | z==pts(2,k) );
        end
        [~,Neighbours] = find(neigh_tf);
        Index = Neighbours;
        Connected(Index)=Mark;
        
    end
    Mark = Mark + Difference;
end

% Compress label values:
v = unique(Connected);
for j = 1:numel(v)
    Connected( Connected==v(j) ) = j;
end

labels = Connected;

%{
clrs = lines(j);
figure,
for j = 1:numel(v)
    plot(zin(:,Connected==j),'Color',clrs(j,:)), 
    hold on, 
end
%}

end


% ------------------------------------------------------------------------
function labels = vdlabel_two_pass(z)
error('not working')
fprintf(2,'See http://en.wikipedia.org/wiki/Connected-component_labeling\n');
%fprintf(2,'Code not yet working properly...\n')
n = size(z,2);
labels = zeros(1,n);
linked = {};
NextLabel = 1;

%labels(1) = 1;

% First pass:
for j = 1:n
    
    % If already labelled, skip,
    if labels(j) ~= 0
        continue
    end
    
    pts = z(:,j);
    z(:,j) = NaN;
    
    % Now find all segments connected to point in top row:
    [~,precol] = find(z==pts(1));
    
    % And all segments connected to point in bottom row:
    [~,postcol] = find(z==pts(2));
    
    %isempty(precol) && isempty(postcol)
    
    % Join:
    neighbourCols = [precol; postcol];
    
    % Check for labelled neighbours 
    neighs = neighbourCols(labels(neighbourCols) > 0);
    %keyboard
    if isempty(neighs)
        % Does not have labelled neighbours:
%         inset = find(cellfun(@(x)any(x==NextLabel),linked));
%         if isempty(inset)
%             linked{end+1} = NextLabel;
%         else
%             keyboard
%             linked{NextLabel} = linked{inset};%
%         end

        %linked{NextLabel} = NextLabel;
        labels(j) = NextLabel;
        labels(neighbourCols) = NextLabel;
        NextLabel = NextLabel+1;
    else
        % Does have labelled neighbours:
        %linked{end+1,1} = [(labels(neighs)) NextLabel];
        linked{end+1,1} = (labels(neighs));
        
        L = labels(neighs);
        labels(j) = min(L);
        %min(labels(lab_neighs))
        keyboard
        for k = 1:numel(L)
            overlap = arrayfun(@(v)any(v==linked{k}),L)
            if any( overlap )
                keyboard
            end
            %arrayfun(@(v)any(v==linked{k}),L)
            %linked{k} = Union(linked,linked{k},L);
        end
    end
end

% Second pass:
lablels_1 = labels;
for j = 1:numel(linked)
    
    fset = find(cellfun(@(x)any(x==j),linked));
    labels(j) = fset(1);
end

keyboard

end %vdlabel_2pass(z)

% ------------------------------------------------------------------------
function L = vdlabel_by_tracing(z)

% Initialise label vector:
L = zeros(1,size(z,2));
l = 1;

j = 0;
while 1
    c = find(L==0,1,'first');
    
    % Exit condition: all zeros are converted to numbers:
    if isempty(c)
        break
    end
    
    % Call the tracing function to trace out with current label value:
    [z,L] = trace_and_label(z,L,c,l);
    
    l = l+1;
    
    % Safety check - remove when working
    j = j+1;
    if j > 5000
        fprintf(2,'infinite loop...\n');
        keyboard
    end

end

end

% ------------------------------------------------------------------------
function [z,L] = trace_and_label(z,L,col,val)
% recursive function
%     z: complex voronoi diagram
%     L: label vector
%   col: column for seeding the start
%   val: label id to insert into label vector L

% Label this point with the current label id:
L(col) = val;

% Store the line segment & remove it from the diagram:
pts = z(:,col);
z(:,col) = NaN;

% Now find all points connected to point in top row:
%[~,precol] = find(z==pts(1));

% And all points connected to point in bottom row:
%[~,postcol] = find(z==pts(2));

% Join:
%conn = [precol; postcol];

[~,conn] = find( z==pts(1) | z==pts(2) );

% Keep only valid points for tracing out:
conn = conn(any(~isnan(z(:,conn))));

% Recurse:
for j = 1:numel(conn)
    % Trace out branches from all points:
    [z,L] = trace_and_label(z,L,conn(j),val);
end


end


% ------------------------------------------------------------------------
function [z] = trim_cooked_ends(z)
endids = points_with_occurrence(real(z),imag(z),@(n)n==1);

THRESH = 40;  % [deg]
mask = true(size(z));

for j = 1:numel(endids) % Should be j = [1,2] 
    zj = trace_from_end(z,endids(j));
    dz = diff(zj);
    v1 = dz(1:end-1);
    v2 = dz(2:end);
    ang = angle(dot(v1,v2,2));
    
    % Find breaks by angle threshold
    k = 1 + find(abs(ang) > THRESH*pi/180,1,'first');
    if k < min( numel(zj)/8 )
        % use points 1:k in (x,y) to drop points from (vx,vy)
        for q = 1:k
            [~,c] = find( z == z(q) );
            mask(:,c) = false;
        end
    end
end

% apply the mask:
z = z(mask);


end %trim_crooked_end()


% ------------------------------------------------------------------------
function [zout] = trace_from_end(zin,id)
% Single-direction trace from endpoint at specified index
j = 0;
[r,c] = ind2sub(size(zin),id);

otherrow = @(rowid)(~(rowid-1))+1;

% Initialise vectors with first endpoint:
zout = zin(r,c);

nmax = numel(zin);

while ~isempty(zin)
    
    % Find second point on line:
    z2 = zin(otherrow(r),c);
    
    % Store:
    zout(end+1,:) = z2;
    
    % Now drop this line:
    zin(:,c) = [];
    
    % This point is also the first point on another line, find it:
    [r,c] = find(zin==z2);
    
    % Done.  
    
    % Just an error catch:
    j = j+1;
    if j > nmax
        fprintf(2,'Cicular problem...\n');
        keyboard
    end     
end

end

% ------------------------------------------------------------------------
function [vx,vy] = erode_voronoi_branches(vx,vy)
% Inputs & outputs are in the form that voronoi produces: 2-by-N matrices

% Trace brances back to their join:
tic
branches = trace_branches(vx,vy);
toc


% % Convert branches to conglomorate voronoi format:
% siz = cell2mat(cellfun(@size,branches,'uniformoutput',false));
% nlines = sum(siz(:,1)) - numel(siz(:,1));
% bvx = NaN(2,nlines);
% bvy = bvx;
% c = 0;
% for j = 1:numel(branches)
%     for k = 1:size(branches{j},1)-1
%         bvx(1:2,c+k) = branches{j}(k:k+1,1);
%         bvy(1:2,c+k) = branches{j}(k:k+1,2);
%     end
%     c = c+k;
% end

% Sort in order of shortest brances by geometric length
len = @(b)sum(sqrt(sum(diff(b).^2,2)));
[~,idx] = sort(cellfun(len,branches));
branches = branches(idx);

% Branches are defined in junction-last order
% Progressively drop branch segments. Drop if root node still has more than
% 2 branches
for j = 1:numel(branches)
    ns = size(branches{j},1);
    bx = branches{j}(:,1);
    by = branches{j}(:,2);
    % Test conditions for dropping a branch:
    r1 = sum(any(vx==bx(1) & vy==by(1)));       % repeats of point1
    r2 = sum(any(vx==bx(end) & vy==by(end)));   % repeats of point2
    if (r1<=1 && r2<=1) || ... % both end points are free (disjoint branch or isolated points)
        (r2 >= 3) % There is still a fork here
        for k = 1:ns-1
            col = all( ((vx==bx(k)) | (vx==bx(k+1))) & ((vy==by(k)) | (vy==by(k+1))) );
            vx(:,col) = [];
            vy(:,col) = [];
        end
    else
        % ignore normally joined points, and other problems
    end
    %figure(98), plot(vx,vy,'k'), hold on, plot(bx,by,'r*'), hold off
    %pause(0.1)
end


% Check that there are no points with occurrence 3. If there are, we need
% to recurse this function....

% while 1
%     p3 = points_with_occurrence(vx,vy,@(n)n>=3);
%     if isempty(p3)
%         break
%     end
% %     
% %     % Initialize
% %     vxk = vx;
% %     vyk = vy;
% %     
% %     % use p3(1) and find its duplicates:
% %     tf = vx == vx(p3(1)) & vy == vy(p3(1));
% %     for j = 1:
% %     
% %     keyboard
% %     
% %     [xk,yk] = trace_from_end(vx,vy,0)
%     
%     fprintf(2,'Unhandled problem...\n');
%     keyboard
% end

%{
figure, plot(vx(1,:),vy(2,:),'k.','MarkerSize',7), grid on, axis equal, hold on
plot(vx(1,:),vy(1,:),'r*')
%}
end


% ------------------------------------------------------------------------
function branches = trace_branches(x,y)
% Trace branches from their end points back to their junction with the
% trunk;
% Branches are defined with the tip first, and the junction last.

% Find the ends of branches, defined by points with 1 occurrence:
[ie,xe,ye] = points_with_occurrence(x,y,@(n)n==1);
siz = size(x); 
branches = cell(numel(ie),1);
otherrow = @(rowid)(~(rowid-1))+1;
for j = 1:numel(ie)
    b = []; % columns of [x y]
    xj = x;
    yj = y;
    
    % Find first point
    [r,c] = ind2sub(siz,ie(j));
    
    b(1,:) = [ xj(r,c) yj(r,c) ];
    b(2,:) = [ xj(otherrow(r),c), yj(otherrow(r),c) ];
    
    % While next point occurs twice, get next connected point, until it is
    % either 3-connected, or 1-connected, then break
    while sum(xj(:) == b(end,1))==2 && sum(yj(:) == b(end,2))==2
        % Find the two items that are the same as b(end):
        [rset,cset] = find(xj == b(end,1) & yj == b(end,2));
         
        % Now get the line data out of the new column:
        c = cset(c~=cset);
        
        % But we need to find which row to get:
        r = rset( b(end,1) ~= x(rset,c) );
        if isempty(r)   % this happens when x-values for start & end of segment are the same
            r = 1;      % it doesn't matter which one we choose then...
        end
        
        % And append onto b:
        b(end+1,:) = [x(r,c) y(r,c)];
        
        % And drop from xy
        x(:,c) = NaN;
        y(:,c) = NaN;
        
    end
    branches{j} = b;    % store for output
end

%{
figure
axes;
hold on
for j = 1:numel(branches)
    b = branches{j};
    plot(b(:,1),b(:,2),'r-*')
end
axis equal, grid on
keyboard
%}

end %trace_branches(x,y)

% ------------------------------------------------------------------------
function [index,x,y] = points_with_occurrence(x,y,condfun)
% Find points in (x,y) that occur n times according to condfun
%
% Points are repeated if their x and y values are repeated
[~,ix] = occurrencesin(x,condfun);
[~,iy] = occurrencesin(y,condfun);

% A point is repeated only if its x and y values are both repeated:
[index,ia,ib] = intersect(ix,iy);
x = x(ia);
y = y(ib);


end %repeated_points

% ------------------------------------------------------------------------
function [v,index] = occurrencesin(x,condfun,tol)
% Find values in x that occur n times according to the anonymous function
% condition condfun
%
% Usage:
%   occurrencesin(x,@(n)n==1)  
%   occurrencesin(x,@(n)n>=3)
if nargin < 3 || tol == 0
    xr = x;
else
    xr = round(x/tol)*tol;
end
xr = xr(:);
[nr, bin] = histc(xr, unique(xr));
multiple = find(condfun(nr));
index    = find(ismember(bin, multiple));
v = x(index);

end %occurrencesin()

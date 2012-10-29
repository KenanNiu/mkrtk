% ------------------------------------------------------------------------
function X = split_cloud_to_traces(xyz)
% Tricky function.  It takes set of xyz points (N-by-3) and divides them up
% into slices.  The means of doing this is optimising the parameters of a
% plane so it fits the most number of points on it.  All those points
% should then belong to the most populated trace.  From there, we simply
% modify the plane origin (position), and not its direction vectors
% (orientation) and repeat until all points are accounted for.
load temp.mat

r  = size(xyz,1); % # rows
o  = mean(xyz);
% Basic vector projection is like this:
% project = @(a,b) dot(a,b,2)*b/norm(b);

% In matrix form, it looks like this:
% bsxfun(@rdivide, bsxfun(@times,dot(a,b,2),b), norm(b) )

% But we also have to define @norm in the conventional way:
norm = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);

while 1

    % Select some indices to define a random plane:
    idx = select_indices('random');
    
    % Use these points to create a normal:
    p1 = xyz(idx(1),:);
    p2 = xyz(idx(2),:);
    p3 = xyz(idx(3),:);
    
    u = p3 - p2; u = u./norm(u);
    v = p1 - p2; v = v./norm(v);
    
    if 1 - dot(u,v) <= 10*eps
        % Skip this ambiguous combination
        fprintf('skipping!\n')
        continue
    end
    n = cross(u,v); n = n./norm(n);
        
    % Project all points onto the plane normal, in matrix form:
    xp = project_points_onto_line(xyz,repmat(n*200,r,1));
    
    % Then turn these into signed distances:
    d = linePosition3d(xp, [o n]);
    
    % Calculate the spread of the clusters:
    %[spread,k] = mean_cluster_spread(d);
    %fprintf('spread = %f\n',spread)
    
    %k = autocluster1d(d,.1);
    %fprintf('k = %d\n',k);
    %if k < 7
    %    keyboard
    %end
    
    
    % Number of slices estimated using this plane orientation:
    tol = 1e-10;
    [c,ia,ib] = unique( round(d/tol)*tol );
    
    % Now each element in c represents a possible set of points which lie
    % on a plane.  For it to be a valid set, there must be at least 3
    % points which define the plane, so 
    cl_ids = unique(ib);
    ok = true;
    j = 0;
    while ok && j < numel(cl_ids)
        j = j+1;
        sel = ib==cl_ids(j);
        ok = sum(sel) >= 3;
    end
    % If all the unique sets contain at least 3 points, we have a possible
    % solution.  It isn't a formally watertight solution, but it'll do for
    % now.
    if j == numel(cl_ids)
        break
    end
        
end %while

% Now we can divide the points up into their groups:
ns = numel(cl_ids);
X = cell(ns,1);
for j = 1:ns
    X{j} = xyz( ib == j , : );
end 

% Done!

    % ---------------------------------
    function inds = select_indices(type)
        switch type
            case 'random'
                % randomly select 3 unique indices:
                inds = round( rand(10,1) * (r-1)) + 1;
                inds = unique(inds);
                inds = inds(1:3);
            case 'predict1'
            case 'predict2'
        end
        
    end

    % ----------------------------------
    function xp = project_points_onto_line(pts,vec)
        xp = bsxfun(@rdivide, bsxfun(@times,dot(pts,vec,2),vec), norm(vec) );
    end
        
%{
figure, plotcloud(xyz,'r*'); hold on, 
quiver3(o(1),o(2),o(3),n(1),n(2),n(3),50), 
grid on, axis equal
o = xyz(idx(1),:);
u = u;
v = cross(u,n); v = v/norm(v);
hp = drawPlane3d([o u v]);
set(hp,'FaceAlpha',0.6)
%hold on, plotcloud(xp,'g*')
keyboard
%}


end %split_cloud_to_traces()

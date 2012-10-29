% ------------------------------------------------------------------------
function X = split_cloud_to_traces(xyz)
% Tricky function.  It takes set of xyz points (N-by-3) and divides them up
% into slices.  The means of doing this is optimising the parameters of a
% plane so it fits the most number of points on it.  All those points
% should then belong to the most populated trace.  From there, we simply
% modify the plane origin (position), and not its direction vectors
% (orientation) and repeat until all points are accounted for.
load temp.mat

% Step 1: Defining the plane orientation
o = mean(xyz,1); u = [1 0 0]; v = [0 1 0];
p = [o u v];
n = cross(u,v);
x0 = [o n];     % [origin, normal]

opts = optimset('fminsearch');
%[x,fval,exitflag,output] = fminsearch(@fun0,x0,opts)
%[x,fval,exitflag,output] = fminunc(@fun0,x0);

mx = max(xyz);
mn = min(xyz);

lb = [mn, 0 0 0];
ub = [mx, 1 1 1];
A = [1 0 0 0 0 0    % ox <= max(x)
    -1 0 0 0 0 0    
    0  1 0 0 0 0    % oy <= max(y)
    0 -1 0 0 0 0
    0 0  1 0 0 0    % oz <= max(z)
    0 0 -1 0 0 0
    0 0 0  1 0 0    % x <= 1
    0 0 0 -1 0 0    % x >= -1
    0 0 0 0  1 0
    0 0 0 0 -1 0
    0 0 0 0 0  1
    0 0 0 0 0 -1
    ];
b = [mx(1)
    mn(1)
    mx(2)
    mn(2)
    mx(3)
    mn(3)
    1
    1
    1
    1
    1
    1]

[x,fval,exitflag,output] = fmincon(@fun0,x0,A,b);

keyboard
% Step 2: Iteratively select the next largest group of points to define a
% trace
%   Actually, fminsearch tends to get stuck on local minimia, so we may not
%   find the slice with the most points in the first instance, but that's
%   ok.

% Step 3: Sort in order of increasing value along the plane normal.

    function n = fun0(x)
        % Minimisation function
        oi = x(1:3);    % Plane origin
        ni = x(4:6);    % Plane normal
        r  = size(xyz,1); % # rows
        
        % Now calculate the normal distance of all points to the plane
        
        
        % Basic vector projection is like this:
        % project = @(a,b) dot(a,b,2)*b/norm(b);
        
        % In matrix form, it looks like this:
        % bsxfun(@rdivide, bsxfun(@times,dot(a,b,2),b), norm(b) )
        
        % But we also have to define @norm in the normal way:
        norm = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);
        
        % Then project all points (a) onto the plane normal (b)
        av = xyz;
        bv = repmat(ni,r,1);
        
        xp = bsxfun(@rdivide, bsxfun(@times,dot(av,bv,2),bv), norm(bv) );
        
        % Now count the number of points that have zero distance to the
        % plane origin (ie, they lie on the plane)
        d3 = xp - repmat(oi,r,1);
        if any(d3) < 0
            keyboard
        end
        n = sum(norm(d3) <= 0.0001);
                
        % And because we're doing minimisation, we negate:
        n = r-n;
        fprintf(' %d',x);
        fprintf('\n');
        
    end %fun


end
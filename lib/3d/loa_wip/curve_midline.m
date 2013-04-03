function M = curve_midline(X)
% X is an N-by-1 cell array of M-by-3 sets of points which represent
% different curves

if iscell(X)
    % Recursion:
    M = cell(size(X));
    for j = 1:numel(X)
%        print_status()
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
    M2 = smooth_midline2d(X2);
    %M2 = skel_midline2d(X2);
                
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

%     % ---------------------
%     function print_status()
%         ds = dbstack;
%         nt = sum(strcmp({ds.name},mfilename))-1;
%         indent = repmat('  ',1,nt);
%         fprintf('%sCalculating midline %d of %d\n',indent,j,numel(X));
%     end

end %curve_midline()


% ------------------------------------------------------------------------
function m = smooth_midline2d(xy)
% first we create an representative image which describes the polygon:

[B,T,S] = poly2bw(xy);

% Internal distance map:
%D = bwdist(~B);

% Smoothing routine: slightly simpler & quicker than an active contour:
inds = find(B);
[yi,xi] = ind2sub(size(B),inds);

% Two-step process quicker than using only 'lowess', & smoother than using
% only 'moving':
span = round(length(B)/5);          % ~20% span
span = span + mod(span-1,2);        % Span must be odd for 'moving'
yi = smooth(xi,yi,span,'moving');   % (1) Moving average filter (fast)
[xu,indu] = unique(xi);             %\_ reduce data down to unique line
yu = yi(indu);                      %/
yu = smooth(xu,yu,0.05,'lowess');   % (2) Smooth curve 

m = [xu yu]*S + repmat(T(1,:),numel(xu),1);     % Re-scale back to original

end %smooth_midline2d()


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

end %skel_midline2d()


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

end %poly3bw()

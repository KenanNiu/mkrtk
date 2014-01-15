% function closest_point_curve_line
% %
% %
% 
% [x,y,z] = create_circular_path(10);
% curve = [x(:), y(:), z(:)];
% 
% 
% % -------- Plotting nonsense ---------
% figure, plot3(x,y,z,'LineStyle','-','Marker','*');
% axis equal
% grid on
% hold on
% title(gca,'Please click to create a line parallel to the view')
% set(gcf,'WindowButtonDownFcn','uiresume')
% % -------------------------------------
% 
% while 1
%     uiwait(gcf)
%     cp = get(gca,'CurrentPoint');
%     
%     edge = [cp(1,:), cp(2,:)];
%     
%     
%     plot3( cp(:,1),cp(:,2), cp(:,3), 'r*-')
%     
%     
%     line = edge2line(edge);
%     
%     [pt_curve, pt_line] = core(curve, line)
%     
%     plot3(pt_curve(1),pt_curve(2),pt_curve(3),'m*')
% end
% 
% function [x,y,z] = create_circular_path(radius)
% theta = linspace(-3*pi/4,pi/4,10);
% x = 0.*theta;
% y = radius.*sin(theta);
% z = radius.*cos(theta);
% 
% 
% % =====================================================

function [pt_curve,pt_line] = core(curve,line)
% CURVE is N-by-3 polyline
% LINE is 2-by-3 line segment: [x1 y1 z1; x2 y2 z2]

% Convert polyline into list of segments:
seg = [curve(1:end-1,:) curve(2:end,:)];

% Unit direction vectors of the segments:
v = seg(:,4:6) - seg(:,1:3);
vunit = normalizeVector3d(seg(:,4:6) - seg(:,1:3));

% Describe infinite rays of each segment in the parameterized form:
%   p0 + t*v
p0 = seg(:,1:3);

seg_rays = [p0 vunit];

[d,pts_line,pts_seg_rays] = mutual_perpendicular(line,seg_rays);

%% Parametrized distance of the closest point on each segment along their
% respective segments:
t = linePosition3d(pts_seg_rays,edge2line(seg_rays));

% Constrain projections to withing bounds of segments (ie, 0 < t < 1)
t( t < 0 ) = 0;
t( t > 1 ) = 1;

% Coordinates of points on the segments:
pt_seg_constrained = bsxfun(@plus, p0, [t.*v(:,1) t.*v(:,2) t.*v(:,3)]);

% Distance between points of each solution:
d = sqrt( sum(pts_line - pt_seg_constrained,2).^2 );

[dmin,id] = min(d);

pt_curve = pt_seg_constrained(id,:);
pt_line  = pts_line(id,:);

% %%
% plot3(pts_line(:,1),pts_line(:,2),pts_line(:,3),'m*')
% %%
% plot3(pts_seg_rays(:,1),pts_seg_rays(:,2),pts_seg_rays(:,3),'g*')
% %%
% plot3(pt_seg_constrained(:,1),pt_seg_constrained(:,2),pt_seg_constrained(:,3),'g*')


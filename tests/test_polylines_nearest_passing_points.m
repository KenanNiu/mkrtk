function test_polylines_nearest_passing_points

% Curve 1
radius = 10;
theta = linspace(-3*pi/4,pi/4,10);
x1 = 0.*theta;
y1 = radius.*sin(theta);
z1 = radius.*cos(theta);

% x1 = x1(4:4);
% y1 = y1(4:4);
% z1 = z1(4:4);

% Curve 2
radius = 50;
theta = linspace(3*pi/4, pi/2,20);
x2 = radius.*sin(theta) - radius;
y2 = radius.*cos(theta);
z2 = 0.*theta + 2;


% x2 = x2(15:16);
% y2 = y2(15:16);
% z2 = z2(15:16);


figure
h1 = plot3(x1,y1,z1,'b-','Marker','*');
hold on
h2 = plot3(x2,y2,z2,'m-','Marker','*');
grid on
axis equal
axis tight


curve1 = [x1(:) y1(:) z1(:)];
curve2 = [x2(:) y2(:) z2(:)];


[p1,p2] = polylines_nearest_passing_points(curve1,curve2);


% 
% edge1 = reshape(curve1',1,6);
% edge2 = reshape(curve2',1,6);
% 
% [p1,p2] = edge_edge_closest_point(edge1,edge2)


plot3(p1(1),p1(2),p1(3),'bo')
plot3(p1(1),p1(2),p1(3),'bx')

plot3(p2(1),p2(2),p2(3),'mo')
plot3(p2(1),p2(2),p2(3),'mx')
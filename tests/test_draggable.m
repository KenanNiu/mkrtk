function test_draggable
close all

w = 10;

[mpath.x,mpath.y,mpath.z] = create_circular_path(w);

figure
curve = plot3(mpath.x,mpath.y,mpath.z+2,'.-');
set(curve,'HitTest','off');
axis equal
grid on
hold on

id = round(2/3*numel(mpath.x));
pt.x = mpath.x(id);
pt.y = mpath.y(id);
pt.z = mpath.z(id);

ln.x = (-5:5) + pt.x + 0;
ln.y = repmat(pt.y,size(ln.x)) + 2 ;
ln.z = repmat(pt.z,size(ln.x)) - 3 ;

point = plot3(pt.x,pt.y,pt.z,'Marker','o','MarkerFaceColor','g','MarkerSize',12);

line = plot3(ln.x,ln.y,ln.z,'Marker','*','Linestyle','-','Color','r');


%draggable(point)
draggable(point,'ConstrainTo',curve)
%draggable(point,'ConstrainTo',curve,'ButtonUpFcn',@my_button_up_fcn)
%draggable(point,'ConstrainTo',curve,'EndFcn',@()my_end_fcn(gca))

%draggable(line)
%draggable(line,'ConstrainTo',curve)
%draggable(line,'AllowRotate',true)
draggable(line,'ConstrainTo',curve,'AllowRotate',true)

axis auto
axis tight


%------------ Interactive plane ----------%
p0  = [0 0 0];
nml = [1 0 0];
p = createPlane(p0,nml);

hp = drawPlane3d(p);
set(hp,'FaceColor',[0.5 0.5 0.5])

draggable(hp,'AllowRotate',true)
%draggable(hp,'AllowRotate',true,'ConstrainTo',curve)

%keyboard


function [x,y,z] = create_circular_path(radius)
theta = linspace(-3*pi/4,pi/4,10);
x = 0.*theta;
y = radius.*sin(theta);
z = radius.*cos(theta);

function my_button_up_fcn(hobj)
title(ancestor(hobj,'axes'),num2str(rand(1)))


function my_end_fcn(ax)
title(ax,'FUNCTION DISABLED!!!!')

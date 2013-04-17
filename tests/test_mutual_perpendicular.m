clear variables

% Test mutual_perpendicular against momentArm


p1 = [...
    0 0 0;
    1 1 1;
    -2 3 4  
    2 4 6]; 
p2 = [...
    5 -0 6
    -4 3 -5
    1 1 1
    -4 5 3]; %this
p3 = [...
    6 5 7
    -9 -4 0
    3 -4 6
    4 -4 -5];
p4 = [...
    5 1 4
    3 3 -4
    9 0 1   % this
    -7 5 -9];

edge1 = [p1 p2];
edge2 = [p3 p4];

line1 = edge2line(edge1);
line2 = edge2line(edge2);


[dmp,P1,P2,sc,tc] = mutual_perpendicular(line1,line2);



figure, 
plot3(edge1(:,[1 4])',edge1(:,[2 5])',edge1(:,[3 6])','color',[0 .6 0])
hold on
plot3(edge2(:,[1 4])',edge2(:,[2 5])',edge2(:,[3 6])','color','r')

plot3([P1(:,1) P2(:,1)]',[P1(:,2) P2(:,2)]',[P1(:,3) P2(:,3)]','Color','k')

grid on, 
axis equal

for j = 1:size(p1,1)
    xyz1 = line1(j,1:3);
    xyz2 = line2(j,1:3);
    text(xyz1(1),xyz1(2),xyz1(3),sprintf('(%d)-1',j))
    text(xyz2(1),xyz2(2),xyz2(3),sprintf('(%d)-2',j))
end
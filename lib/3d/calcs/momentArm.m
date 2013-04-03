function [d,P1,P2] = momentArm(p1,p2,p3,p4)
%MOMENTARM Calculate moment arm from two lines, each defined by two points.
%
% See also distanceLines3d.m in geom3d toolbox


%==== METHOD 1 ====%
% See :
%   http://geomalgorithms.com/a07-_distance.html
%   http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm

TOL = 1e-7;
norm = @(x)sqrt(dot(x,x));

u = p2-p1;
v = p4-p3;
w = p1-p3;

a = dot(u,u);
b = dot(u,v);
c = dot(v,v);
d = dot(u,w);
e = dot(v,w);
D = a*c-b*b;

if D < TOL
    sc = 0;
    tc = max( d/b, e/c );
else
    sc = (b*e - c*d) / D;
    tc = (a*e - b*d) / D;
end

d3 = w + (sc*u) - (tc*v);

d = norm(d3);

if nargout > 1
    % Intersection is defined at these two points:
    
    % Q(tc) = Q0 + tc*(Q1-Q0)
    Qtc = p3 + tc*v;
    % P(sc) = P0 + sc*(P1-P0)
    Psc = p1 + sc*u;
    
    % Outputs:
    P1 = Qtc;
    P2 = Psc;
end
%}

%==== METHOD 2 ====%
%{
% http://www.coventry.ac.uk/ec//jtm/slides/8/sld8p5.pdf

% The mutual perpendicular vector is the cross product:
v3 = cross(v1,v2)

% The shortest distance:
d = abs( dot((p2-p1), v3) / norm(v3) )
%}




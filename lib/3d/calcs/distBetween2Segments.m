function [distance,point1,point2] = distBetween2Segments(edges1,edges2)
%DISTBETWEEN2SEGMENTS computes the minimum distance between two line
% segments and the corresponding endpoints of the shortest joining segment.
%
% Adapted from  from Dan Sunday's Geometry Algorithms originally written  
% in C++:
% http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment
%
% Code has been heavily modified to support full vectorisation.
% 
% Inputs are two sets of line segments, where each one is a row
% describing the start and end points of the segment: [x1 y1 z1 x2 y2 z2] 
%
% If EDGES1 and EDGES2 have the same number of rows, the closest points are
% calculated for each edge pair: edges1(i,:) , edges2(i,:)
% Otherwise, either EDGES1 or EDGES2 should have only one row and the
% nearest points are calculated for each of the other edges, such that P1
% and P2 are both N-by-3.
%
% Joshua Martin, 12-Apr-2013

[edges1,edges2] = conformsize(edges1,edges2);

p1 = edges1(:,1:3);
p2 = edges1(:,4:6);
p3 = edges2(:,1:3);
p4 = edges2(:,4:6);

u = p2 - p1;
v = p4 - p3;
w = p1 - p3;

a = dot(u,u,2);
b = dot(u,v,2);
c = dot(v,v,2);
d = dot(u,w,2);
e = dot(v,w,2);
D = a.*c - b.*b;

sD = D;
tD = D;
sN = NaN(size(sD));
tN = NaN(size(sN));

SMALL_NUM = 0.00000001;

%=== compute the line parameters of the two closest points === %

% -- Deal with parallel lines
ip = D < SMALL_NUM; % these lines are almost parallel
np = ~ip;
sN(ip) = 0.0;       % force using point P0 on segment S1
sD(ip) = 1.0;       % to prevent possible division by 0.0 later
tN(ip) = e(ip);
tD(ip) = c(ip);

% -- Deal with non-parallel lines:

% get the closest points on the infinite lines
sN(np) = (b(np).*e(np) - c(np).*d(np));
tN(np) = (a(np).*e(np) - b(np).*d(np));

    idx = np & (sN < 0);    % sc < 0 => the s=0 edge is visible
    sN(idx) = 0;
    tN(idx) = e(idx);
    tD(idx) = c(idx);

    idx = np & (sN > sD);   % sc > 1 => the s=1 edge is visible
    sN(idx) = sD(idx);
    tN(idx) = e(idx) + b(idx);
    tD(idx) = c(idx);

% ... t=0 visible edges ... % tc < 0 => the t=0 edge is visible
vis = tN < 0;
tN(vis) = 0.0;

    % recompute sc for these edges
    idx1 = vis & (-d < 0);
    sN(idx1) = 0;

    idx2 = vis & (-d > a);
    sN(idx2) = sD(idx2);

    idx3 = vis & ~(idx1 | idx2);

    sN(idx3) = -d(idx3);
    sD(idx3) = a(idx3);


% ... t=1 visible edges ... % tc > 1 => the t=1 edge is visible
vis = (tN > tD);
tN(vis) = tD(vis);

    % recompute sc for these edges
    idx1 = vis & ( (-d+b) < 0 );
    sN(idx1) = 0;

    idx2 = vis & ( (-d+b) > a );
    sN(idx2) = sD(idx2);

    idx3 = vis & ~(idx1 | idx2);
    sN(idx3) = -d(idx3) + b(idx3);
    sD(idx3) = a(idx3);


% ======= finally do the division to get sc and tc ======= %
tf = abs(sN) < SMALL_NUM;
sc = zeros(size(tf));   % Initialize
%sc(tf) = 0;            % Now redundant
sc(~tf) = sN(~tf)./sD(~tf);

tf = abs(tN) < SMALL_NUM;
tc = zeros(size(tf));       % Initialize
%tc(tf) = 0;                % Now redundant
tc(~tf) = tN(~tf)./tD(~tf);


sc_x_u = bsxfun(@times,sc,u);
tc_x_v = bsxfun(@times,tc,v);

dP = w + sc_x_u - tc_x_v;


distance = sqrt(sum(dP.^2,2));

point1 = edges1(:,1:3) + sc_x_u;
point2 = edges2(:,1:3) + tc_x_v;

end



% ---------------------------------------
function [e1,e2] = conformsize(e1,e2)
[nr1,nc1] = size(e1);
[nr2,nc2] = size(e2);

assert( nc1==6 && nc2==6, 'Inputs must have 6 columns')

% Explicitly conform size:
if nr1 == nr2
    % Fine
elseif nr1 == 1
    e1 = repmat(e1,nr2,1);   % conform size of edges1 to edges2
elseif nr2 == 1
    e2 = repmat(e2,nr1,1);   % conform size of edges2 to edges1
else
    error('incorrectly sized input arguments')
end

end
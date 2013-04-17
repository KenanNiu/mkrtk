function [d,p1,p2] = mutual_perpendicular(line1,line2)
%MUTUAL_PERPENDICULAR Calculate shortest joining line segment (mutual
% perpendicular segment) between two lines or sets of lines.
%
% The inputs LINE1 & LINE2 can be 1-by-6 or N-by-6 lists of line 
% definitions: [x0 y0 z0 u v w]
%
% If one input is 1-by-6 and the other is N-by-6, then d will be N-by-1.
% If both inputs have more than one row, then the number of rows must be
% the same.

[line1,line2] = parse_inputs(line1,line2);

% Method, see :
%   http://geomalgorithms.com/a07-_distance.html
%   http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm

%norm = @(x)sqrt(dot(x,x));

p0 = line1(:,1:3);
q0 = line2(:,1:3);

u = line1(:,4:6);
v = line2(:,4:6);
w = line1(:,1:3)-line2(:,1:3);

a = dot(u,u,2);
b = dot(u,v,2);
c = dot(v,v,2);
d = dot(u,w,2);
e = dot(v,w,2);
D = a.*c-b.*b;


% Deal with parallel lines:
TOL = 10*eps;

% --- Non-vectorised form --- %

% if D < TOL
%     sc = 0;
%     tc = max( d./b, e./c );
% else
%     sc = (b.*e - c.*d) ./ D;
%     tc = (a.*e - b.*d) ./ D;
% end

% --- The above in vectorised form: --- %

ip = (D < TOL); % is parallel
np = ~ip;       % not parallel

sc(ip) = 0;
maxof = max( d./b, e./c );
tc(ip) = maxof(ip);

sc(np) = (b(np).*e(np) - c(np).*d(np)) ./ D(np);
tc(np) = (a(np).*e(np) - b(np).*d(np)) ./ D(np);

% ---   --- %

sc_x_u = bsxfun(@times,u,sc);
tc_x_v = bsxfun(@times,v,tc);

d3 = w + sc_x_u - tc_x_v;%w + (sc.*u) - (tc.*v);

d = sqrt(sum(d3.^2,2));

if nargout > 1
    % Intersection is defined at these two points:
    
    % P(sc) = P0 + sc*(P1-P0)
    Psc = p0 + sc_x_u;          %p0 + sc*u;
    % Q(tc) = Q0 + tc*(Q1-Q0)
    Qtc = q0 + tc_x_v;          %q0 + tc*v;
    
    % Outputs:
    p1 = Psc;
    p2 = Qtc;
end

% ------------------------------------------------
function [line1,line2] = parse_inputs(line1,line2)
[nr1,nc1] = size(line1);
[nr2,nc2] = size(line2);

assert( nc1==6 && nc2==6, 'Inputs must have 6 columns')

if all([nr1 nr2] > 1) && ( nr1 ~= nr2)
    error('Inputs must have the same number of rows if numrows > 1')
elseif (nr2 == 1) && (nr1 > 1)
    line2 = repmat(line2,nr1,1);
elseif (nr1 == 1) && (nr2 > 1)
    line1 = repmat(line1,nr2,1);
end

    
    

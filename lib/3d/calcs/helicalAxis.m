function [axs,theta,pt,slide] = helicalAxis(A1,A2,B1,B2)
%HELICALAXIS Calculate the helical axis of object B with respect to A
%
% [axs,ang,pt,slide] = helicalAxis(A1,A2,B1,B2) calculates the
% instantaneous helical axis of object B with respect to object A.  
%
% As object A moves from position A1 to A2, object B also undergoes motion
% from B1 to B2 in the global space. The motion of B wrt A is found by
% removing object A's motion from object B, leaving a relative motion which
% is equivalent to moving from B1 to B1_rel.  The helical axis returned is the
% instantaneous axis of rotation required to rotate object B from B1 to
% B1_rel.
%
% Inputs
%   - Bone A: N-by-3 matrix of points
%   - Bone B: M-by-3 matrix of points


% PROCRUSTES determines the least squares mapping of Y to X and produces
% the rotation matrix and the translation vector.  These define, in order:
%   1) rotation about the world (global) axis system, then
%   2) translation to the new position 
%
% PROCRUSTES solves the equation:
%   Z = b*Y*T + c;
% where
%   Z = X + err
%
% [E,Z,transform] = procrustes(X,Y);

%%
%{
%% Hashizumi
C2 = A2.CS*B2.CS';
C1 = A1.CS*B1.CS';
R = C1'*C2;
%% Sometimes sames as procrustes
Rb = (B1.CS*B2.CS');
Ra = (A1.CS*A2.CS');
R_8 = Ra*Rb;
%% Other times same as procrustes
Rb = (B2.CS*B1.CS');
Ra = (A2.CS*A1.CS');
R_9 = Rb'*Ra';      
%}

%
%% Using CS and calculations at position (1)
%   Position (3) is position (2) corrected for movement in A
RA12 = A2.CS*A1.CS';
RB12 = B2.CS*B1.CS';
RB13 = RB12*RA12';         % Same as legacy
R = RB13;

TA12 = (A2.Origin' - RA12*A1.Origin')';
TA21 = -TA12;
TB13 = RA12'*B2.Origin' + TA21' - RB13*B1.Origin';
d = TB13';                 % Not exactly same as legacy, but very similar
%}
%%
%{
%% Legacy method - procrustes, and calculations at (2):
% Determing transformation to take A1 to A2
[~,~,tformA] = procrustes(A2.xyz, A1.xyz, 'scaling', false, 'reflection', false);
RA = tformA.T';     % Procrustes uses post-multiplication, hence transpose
TA = tformA.c;

% Remapping TA to the size of B1
TA_sizeB=repmat(TA(1,:),[size(B1.xyz,1) 1]);

% Performing this transformation on B1, to take it to Brel 
Brel_xyz = (RA*B1.xyz' + TA_sizeB')';

% Determining Transformation to take Brel to B2 
[~,~,tformB] = procrustes(B2.xyz, Brel_xyz, 'scaling', false, 'reflection', false);

R = tformB.T';      % Procrustes uses post-multiplication, hence transpose
d = tformB.c(1,:);
%}
%%

% R is now used to find the axis and angle of rotation.  There are two
% ways of doing this:
%   1) Extract axis & angle from rotation matrix
%   2) Convert R to quaternion, convert quaternion to axis/angle

% Method 1):
%   a) Use inbuilt matlab function, if the user has it installed
if exist('vrrotmat2vec', 'file') == 2
    v = vrrotmat2vec(R);
end
%theta = v(4);    % Finding the angle of rotation 
%axs = v(1:3);    % Finding unit vector in direction of axis of rotation 
% Same as:
%   b) Springer Handbook of Robotics, Springer, 2008, p.16
%       http://books.google.com.au/books?id=Xpgi5gSuBxsC&pg=PA16
sign = @(x)(x>=0)*2-1;  % same as Matlab's sign(x), but assigns sign(0)==+1;
m = [R(3,2) - R(2,3), R(1,3) - R(3,1), R(2,1) - R(1,2)];
theta = sign(m*d') * abs( acos( (trace(R)-1)/2 ) );
w = m/(2*sin(theta));
rho = (eye(3) - R')*d'/(2*(1-cos(theta))); % <- point on the axis
rho = rho';

% Method 2):
[theta,axs] = quaternion.RotationMatrix(R).angleaxis;
axs = axs(:)';
            
% Rodrigues vector b
b = tan(theta/2)*axs;

% Translation perpendicular to axis
dstar = d + b;

% Slide along axis:
slide = axs*d';   %Equivalent: slide = dot(d,axs); 

% Determining fixed point of the spatial displacement
%   (From: Geometric Design Of Linkages, McCarthy, Springer 2010)
%   (equivalent to Woltring et. al, 1985 ?)
pt_star = (cross(b, (dstar - cross(b,dstar))))/(2*dot(b,b));

% actually because cross(b,d) = cross(b,dstar), we can just do
pt = (cross(b, (d - cross(b,d))))/(2*dot(b,b));


% Carry out checks if matlab function has also been used:
if exist('v','var') % calculated above with vrrotmat2vec
    if ( max( abs(axs - v(1:3)) ) > 1e12 ) ||...
       ( max( abs(axs - w) ) > 1e12 )
        fprintf('Axis solutions differ...\n')
        keyboard
    end
    if max( abs(pt - pt_star) ) > 1e12
        fprintf('Point solutions differ...\n')
        keyboard
    end
end


% ---------- formal component complete -------- %


% The slide vector is now defined by pt and axs.  However the position of
% pt along axs is somewhat arbitrary, so let's slide it along axs to the
% point that is nearest the cetroid of A & B.
% (see linePosition3D in geom3d by David Legland on FEX)
CG = mean([A1.xyz; A2.xyz; B1.xyz; B2.xyz]);

% Vector from line origin to CG:
vcg = bsxfun(@minus, CG, pt);

% Compute position using dot product normalized with norm of line vector:
t_recenter = bsxfun(@rdivide, sum(bsxfun(@times, vcg, axs), 2), sum(axs.^2, 2) );

% New point:
pt = pt + t_recenter*axs;


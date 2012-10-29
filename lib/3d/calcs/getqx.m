function [q,x] = getqx(B,A)
%GETQX Determine the rotation quaternion and translation that maps A to B.
%
% T = GETQX(B,A) first solves the homogeneous equation:
%   B = TA
% which transforms A to B with the homogeneous transformation T.  This uses
% a Statistics Toolbox function, PROCRUSTES to determine the least squares
% mapping of A to B.   
% PROCRUSTES solves the equation:
%   Z = b*Y*T + c;
% where
%   Z = X + err
%
% [d,Z,transform] = procrustes(X,Y);
%
% Using the results in TRANSFORM, we calculate the rotation quaternion, and
% use the translation component directly
%

[~,~,tform] = procrustes(B, A, 'scaling', false, 'reflection', false);  

q = quaternion.RotationMatrix(tform.T);
x = tform.c(1,:);

% Now:
%   B == q.rotate(A) + X


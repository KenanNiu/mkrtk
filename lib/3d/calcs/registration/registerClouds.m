function [q,x] = registerClouds(solver,static0,dynamic,qi,xi)
% Register static cloud with dynamic cloud
%
% registerClouds('solver',s0,dynamic,q,x)

% Make pre-positioned cloud:
s1 = static0.transform(qi,xi);

% Calculate transformation which maps s1->dynamic:
[R,T] = register(dynamic.xyz,s1.xyz,solver);

% Now recover the transformation which maps s0->s2, where s2 is s0 overlaid
% on dynamic:
s2 = s1.transform(R,T);
[q,x] = gettransform(static0,s2,'quat');

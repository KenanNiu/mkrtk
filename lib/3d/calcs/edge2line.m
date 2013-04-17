function line = edge2line(seg)
% Convert list of segements (edges) [x1 y1 z1 x2 y2 z2] to 
% lines: [x0 y0 z0 u v w]

uvw = normalizeVector3d( seg(:,4:6) - seg(:,1:3) );

line = [seg(:,1:3) uvw];
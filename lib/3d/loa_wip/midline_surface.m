function S = midline_surface(M)
% Put a surface through trace midlines

if ~iscell(M)
    error('Curves must come as a set in a cell array: { [n-by-3]; [m-by-3]; ...}')
    
elseif iscell(M) && iscell(M{1})
    
    % Recursion:
    S = cell(size(M));
    for j = 1:numel(M)
        S{j} = midline_surface(M{j});
    end
    
else
    
    % Main function call:
    S = surface_through_points(M);
    
end


end %midline_surface()



% ------------------------------------------------------------------------
function S = surface_through_points(P)
% P is a N-by-1 cell array of M-by-3 point sets
% S is the surface definition with fields:
%   S.x
%   S.y
%   S.z

assert(iscell(P))

%% Conglomerate PCA:
X = cell2mat(P);
coeff = princomp(X);

% Rotate the data so we can use an (x,y) grid with z elevation values.  We
% define a plane using the first principal component and the trace normal.
% We want to define a rotation matrix using the first principal component
% and the trace normal.  However, these are not necessarily orthogonal, so
% we favor the trace normal, and re-calculate an approximate first
% princiapl component which is normal:
r1 = coeff(:,1);
r2 = planeNormal( planeThroughTrace(P{1}) )';
normalize = @(v) v/norm(v);
r3 = normalize( cross(r1,r2) ); % r1 & r2 are not orthogonal
r1 = cross(r2,r3);              % but r2 & r3 are, so use them to get a new r1


R = [r1, r2, r3]';                  % Rotation matrix
xyz = (R*X')';                      % -> Rotate
%t = -mean(xyz1);                    % Translation vector
%xyz = xyz1 + repmat(t,size(xyz1,1),1);    % -> Translate to origin


% Create padded bounding box:
bb = [min(xyz); max(xyz)];
bb = bb(:)';    % [xmin xmax ymin ymax zmin zmax]
pad = min([bb(2)-bb(1), bb(4)-bb(3)])*.15;  % Pad value: 15% of smallest range
bb = bb + [0 0 -pad pad 0 0 ];

% Discretise in x,y with an approximately square grid with a target number of points
ntarget = 1000; % use approximately this many points in the grid (ie, prod([nx,ny]))
dtarget = sqrt( prod( [bb(2)-bb(1), bb(4)-bb(3)] )/ntarget ); 
nx = round((bb(2)-bb(1))/dtarget);
ny = round((bb(4)-bb(3))/dtarget);
nx = max([nx 3]);   % check 
ny = max([ny 3]);   % check
xi = linspace(bb(1),bb(2),nx);
yi = linspace(bb(3),bb(4),ny);

% Fit surface:
zi = gridfit(xyz(:,1),xyz(:,2),xyz(:,3),xi,yi,'smoothness',0.1); 

% Convert to [x,y,z] format:
[r,c] = size(zi);
[xi,yi] = meshgrid(xi,yi);
Xi = [xi(:), yi(:), zi(:)];

% Transform back, in reverse order:
S3 = Xi*(R);

% Then revert to matrix form for the output:
S.x = reshape(S3(:,1),r,c);
S.y = reshape(S3(:,2),r,c);
S.z = reshape(S3(:,3),r,c);

% Ploting
%{
X2 = xyz*(R);
figure, hold on
%plot3(X(:,1),X(:,2),X(:,3),'k+');
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'b+');
grid on, axis equal tight
view(3)

hs = surf(xi,yi,zi);
set(hs,'FaceColor','g','FaceAlpha',0.5)
axis tight

figure,
hs = surf(S.x,S.y,S.z);
set(hs,'FaceColor','g','FaceAlpha',0.5)
hold on, plot3(X(:,1),X(:,2),X(:,3),'bx')
hold on, plot3(X2(:,1),X2(:,2),X2(:,3),'r+')
axis equal tight
%}

end %surface_through_points()

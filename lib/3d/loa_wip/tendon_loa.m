function Loa = tendon_loa(model,studyName,joint)
% Work in progress which calculates the tendon line of action


%Loa = ankle_preparing_for_knee_tendon_loa(model,studyName);
%Loa = ankle_tendon_loa(model,studyName);


%Loa = knee_tendon_loa(model,studyName);
Loa = calc_loa(model,studyName,joint);



global filter_width % defined in smooth_motion>smooth_models
if isempty(filter_width)
    filter_width = 0.99;
end

%% This needs cleaning up and consolidating with smooth_motion()

handles = guidata(gcbf);
theta = joint_dt_proxy(handles.Models);
sfun = @(y)smooth(theta,y,filter_width,'rloess');
for c = 1:size(Loa,2)
    Loa(:,c) = sfun( Loa(:,c) );
end

%%



end

% ------------------------------------------------------------------------
function Loa = knee_tendon_loa(model,studyName)

X = all_tendon_data(model);

handles = guidata(gcbf);
p = get(handles.PhaseSlider,'Value');      % Current phase

% Define an action plane:
msp = get_simple_action_plane();

% Manual edits:
switch studyName
    case 'MGR_015'
        msp(1:3) = msp(1:3) + cross(msp(4:6),msp(7:9))*3;
end

% plot all the data
hf = figure; ax = axes; hold on
for k = 1:numel(X)
    tk = cell2mat(X{k});
    plot3(tk(:,1),tk(:,2),tk(:,3),'k','Marker','.')
end
% Add keypoint of msp:
plot3(msp(1),msp(2),msp(3),'Marker','o','MarkerFaceColor','g');
axis equal tight
grid on, view(3)


% Get midlines of curves:
global M S g_study_name  % <-- temporary measure to speed things up
if isempty(S) || ~isequal(g_study_name,studyName)
    M = curve_midline(X);
    S = midline_surface(M);
    g_study_name = studyName;  % flag for checking if study has changed
else
    disp('Midlines already calculated... referring to those')
end

hs = surf(S{p}.x, S{p}.y, S{p}.z,'Parent',handles.axes1); 
set(hs,'FaceColor',[0.3 0.5 0.6],'FaceAlpha',0.8,'EdgeColor',[0.3 0.5 0.6])
camlight('left'), camlight('right')
lighting gouraud

keyboard

end


% ------------------------------------------------------------------------
function Loa = ankle_preparing_for_knee_tendon_loa(model,studyName)
X = all_tendon_data(model);

% Set upd some graphics stuff:
handles = guidata(gcbf);
%figure, ax = axes; hold on
ax = handles.axes1;
p = get(handles.PhaseSlider,'Value');      % Current phase

fprintf(2,'Just look at current phase; remove this later...\n');
X = X(p);

% Define an action plane, just a simple one for now:
ap = get_simple_action_plane(handles);

% Get midlines of curves:
M = curve_midline(X);

%% Plot all the data
for k = 1:numel(X)
    tk = cell2mat(X{k});
    mk = cell2mat(M{k});
    plot3(ax,tk(:,1),tk(:,2),tk(:,3),'k.','MarkerSize',7)
    plot3(ax,mk(:,1),mk(:,2),mk(:,3),'g.','MarkerSize',7)
end
axis equal tight
grid on, view(3)


%%
S = midline_surface(M);

% Now intersect the surface and the plane
C = intersect_surf_plane(S,ap);

%%
hold on, 
hs = surf(S{1}.x, S{1}.y, S{1}.z); 
set(hs,'FaceColor',[0.3 0.5 0.6],'FaceAlpha',0.5,'EdgeColor',[0.3 0.5 0.6])
plot3(C{1}(:,1),C{1}(:,2),C{1}(:,3),'r-')
%%
fprintf(2,' >> Now we need to put points on the surface to define the tendon\n');
keyboard
% Now we need to do the clipping - either manual/hard coded, or consider
% constrained draggable point/plane.  Use FEX "draggable" for starting
% point:
% http://www.mathworks.com/matlabcentral/fileexchange/4179-draggable

pOffsets = getPlaneOffsets(studyName,'knee');

% Lower truncation plane:
ltp = [ap(1:3) + pOffsets(1)*ap(4:6),  ap(7:9), planeNormal(ap) ];

% Upper Truncation plane:
utp = [ap(1:3) + pOffsets(2)*planeNormal(ltp), ltp(4:end)];


figure, hold on
for k = 1:numel(C)
    plot3(C{k}(:,1),C{k}(:,2),C{k}(:,3),'k-')
end
view(3)
axis auto, axis equal, axis tight

% Show planes:
hltp = drawPlane3d(ltp);
hutp = drawPlane3d(utp);
set([hltp,hutp],'FaceAlpha',0.6)
%%

keyboard

end %knee_tendon_loa()


% ------------------------------------------------------------------------
function actionplane = get_knee_action_plane
% For the knee, doing a conglomerate PCA on the static Femur/tibia actually
% gives us a coronal plane.  So instead let's get all the dynamic clouds of
% all models (including tendon), smash that into one cloud, and do PCA on
% that.
handles = guidata( gcbf );
mdls = handles.Models;
np = max(arrayfun(@(m)numel(m.LoRes),mdls));
for j = numel(mdls):-1:1
    C{j} = (cat(1,mdls(j).LoRes.xyz));
end
XYZ = cat(1,C{:});
PC = princomp(XYZ);

actionplane = [mean(XYZ), PC(:,1)', PC(:,2)'];

end %get_knee_action_plane


% ------------------------------------------------------------------------
function actionplane = get_simple_action_plane()
%{
% Calculate centroid:
XYZ = cell2mat(cat(1,X{:}));
CG = mean(XYZ);

% Conglomerate PCA:
[coeff] = princomp(XYZ);
coeff(:,3) = cross(coeff(:,1),coeff(:,2)); % This ensures a righ-hand axis system, since princomp doesn't care

% Nearest scan plane:
nsp = nearest_scan_plane(X{1},CG);

% Tendon mid-saggital action plane
actionplane = midSaggitalPlane(nsp,CG,coeff);

%}

handles = guidata( gcbf );
fprintf('Warning!! - this will need to be changed to not use HANDLES...\n');
% Get all the models in position 1 & smash into one cloud
p = 1;
n = numel(handles.Models);
C = cell(n,1);
for j = 1:numel(handles.Models)
    mdl = handles.Models(j);
    if isempty(mdl.HiRes)
        xyz = [];
    else
        cld = mdl.HiRes.transform(mdl.q(1),mdl.x(1,:));
        xyz = cld.xyz;
    end
    C{j} = xyz;
end
XYZ = cat(1,C{:});


% Now put plane through all that:
coeff = princomp(XYZ);
coeff(:,3) = cross(coeff(:,1),coeff(:,2)); % This ensures a righ-hand axis system, since princomp doesn't care

cg = mean(XYZ);
actionplane = [cg coeff(:,1)' coeff(:,2)'];


%{
%% Plotting

plot3(CG(1),CG(2),CG(3),'r*')

[~,ind] = max(coeff);
scale = std(XYZ);
scale = scale(ind);
c1 = coeff(:,1);
c2 = coeff(:,2);
c3 = coeff(:,3);
quiver3(CG(1),CG(2),CG(3),c1(1),c1(2),c1(3),scale(1),'r-')
quiver3(CG(1),CG(2),CG(3),c2(1),c2(2),c2(3),scale(2),'r-')
quiver3(CG(1),CG(2),CG(3),c3(1),c3(2),c3(3),scale(3),'r-')

hap = drawPlane3d(actionplane);

%}


end %get_simple_action_plane()


% ------------------------------------------------------------------------
function Loa = calc_loa(model,studyName,joint)

X = all_tendon_data(model);

% Define an action plane:
switch lower(joint)
    case 'ankle'
        msp = get_simple_action_plane();
        
    case 'knee'
        msp = get_knee_action_plane();
        % But hang it close to the centroid of all the tendon data;
        splat = cellfun(@(x)cat(1,x{:}),X,'Uniformoutput',false);
        msp(1:3) = mean( cat(1,splat{:}) );
    otherwise
        error('bad joint name')
end

% plot all the data
hf = figure; ax = axes; hold on
for k = 1:numel(X)
    tk = cell2mat(X{k});
    plot3(tk(:,1),tk(:,2),tk(:,3),'k','Marker','.')
end
% Add keypoint of msp:
plot3(msp(1),msp(2),msp(3),'Marker','o','MarkerFaceColor','g');
axis equal tight
grid on, view(3)


%% Get midlines of curves:
global M S g_study_name  % <-- temporary measure to speed things up
if isempty(S) || ~isequal(g_study_name,studyName)
    M = curve_midline(X);
    S = midline_surface(M);
    g_study_name = studyName;  % flag for checking if study has changed
else
    disp('Midlines already calculated... referring to those')
end

% Now intersect the surface and the plane
C = intersect_surf_plane(S,msp);

% Get main axis of these lines:
Cpc = princomp(cat(1,C{:}));
pc1 = Cpc(:,1)';
pc2 = Cpc(:,2)';

% Get clipping plane offsets
pOffsets = getPlaneOffsets(studyName,joint); 

% Lower truncation plane:
ltp = [msp(1:3) + pOffsets(1)*pc1,  pc2, planeNormal(msp) ];

% Upper Truncation plane:
utp = [msp(1:3) + pOffsets(2)*planeNormal(ltp), ltp(4:end)];

% Draw
hmsp = drawPlane3d(msp);
set(hmsp,'FaceAlpha',0.6,'FaceColor',[0.5 0.5 0.5])

% Show planes:
hltp = drawPlane3d(ltp);
hutp = drawPlane3d(utp);
hmsp = drawPlane3d(msp);
set([hltp,hutp,hmsp],'FaceAlpha',0.6)


%
% Re-plot
figure(hf), clf(hf)
ax = axes; hold on
for k = 1:numel(X)
    
    % Plot tendon curves for this phase:
    tk = cell2mat(X{k});
    %plot3(tk(:,1),tk(:,2),tk(:,3),'-','Marker','.','Markersize',7);
    hd = plot3(tk(:,1),tk(:,2),tk(:,3),'--','Color',[1 1 1]*0.7);
    
    % Plot midlines:
    m = cell2mat(M{k});
    plot3(m(:,1),m(:,2),m(:,3),'g.','Marker','.','Markersize',7)
    
    % Plot midline surface:
    hs = surf(S{k}.x, S{k}.y, S{k}.z); 
    set(hs,'FaceColor',[0.3 0.5 0.6],'FaceAlpha',0.5,'EdgeColor',[0.3 0.5 0.6])

    % Plot intersection curve, C:
    plot3(C{k}(:,1),C{k}(:,2),C{k}(:,3),'r-')
    if k == 1
        axis equal, axis tight,
        grid on, view(3)
    end
    
    title(sprintf('k = %d',k))
    hmsp = drawPlane3d(msp);
    hltp = drawPlane3d(ltp);
    hutp = drawPlane3d(utp);
    set([hmsp,hltp,hutp],'FaceAlpha',0.6)
    
    
    
    pause
    cla(ax)
end
%}

% Now if that all looks good, we can just plot all the intersection curves:
figure(hf), clf(hf)
ax = axes; hold on
for k = 1:numel(X)
    plot3(C{k}(:,1),C{k}(:,2),C{k}(:,3),'r-')
end
grid on, axis equal, axis tight, view(3)


% Plot the truncation planes:
hltp = drawPlane3d(ltp);
hutp = drawPlane3d(utp);
set([hltp,hutp],'FaceAlpha',0.6)
pause

% Truncate:
C = trunc(C,ltp,'above');
C = trunc(C,utp,'below');


%% Define a roughly orhtogonal cutting plane:
mcoeff = princomp(cell2mat(C));
cp_n = mcoeff(:,1);
cp_u = mcoeff(:,2);
cp_v = cross(cp_n,cp_u);
cp_o = mean(cell2mat(C));
cp = [cp_o(:)' cp_u(:)' cp_v(:)'];               % cutting plane


%
% Now for each midline, create a regression curve & intersect that vector
% with the cutting plane to define the line's point and vector:
cla(ax)
Loa = NaN(size(C,1),6);
for k = 1:numel(C)
    ck = C{k};
    lcoeff = princomp(ck);
    loa_v = lcoeff(:,1);  % primary component
    loa_o = mean(ck);        % centroid
    loa = [loa_o(:)', loa_v(:)'];   % Initial definition of the line
    loa(1:3) = intersectLinePlane(loa,cp); % revised anchor point
    
    Loa(k,:) = loa;
    
    plot3(ck(:,1),ck(:,2),ck(:,3),'g-','Marker','.','MarkerSize',7);
    drawLine3d(loa);
    if ~isequal(get(hf,'CurrentKey'),'s')
        pause 
    end
    cla(ax)
end
%}
% Done!

% Loa is defined in rows of (point,vector):
%   [ x0 y0 z0 xhat yhat zhat ]

end %ankle_tendon_loa()


% ------------------------------------------------------------------------
function off = getPlaneOffsets(name,joint)
% Get clipping plane offsets
% Format:
%   { study_name, [upper; lower] }
switch joint
    case 'ankle'
        map = {...
            'Pilot 1', [50; -20]
            'Pilot 2', [55; 0]
            'Pilot 3', [nan; nan]
            'Pilot 4', [55; -5]
            'Pilot 5', [55; 10]
            'Pilot 6', [80; 40]
            'Pilot 7', [70; 30]
            'Pilot 8', [45; 00]
            'MGR_PILOT9', [50; 5]
            'MGR_PILOT10',[65; 25]
            'MGR PILOTT11',[65; 10]
            'MGR PILOTT12',[45; -15]
            'MGR_011',[60 15];
            'MGR_012',[80 45];
            'MGR_013',[70 10];
            'MGR_014',[55 15];
            'MGR_015',[65 15];
            'MGR_017',[70 25];
            'MGR_018',[80 30];
            'MGR_019',[50 20];
            'MGR_028',[60 15];
            'MGR_033',[95 55];
            'Josh_s_stuff_Data_11_July_14',[50 0];
            };
    case 'knee'
        map = {...
            'MGR_011', [40; 20]
            'MGR_012', [25; -15]
            'MGR_013', [10; -15]
            'MGR_014', [45; 20]
            'MGR_015', [40; 10]
            'MGR_017', [30; -10]
            'MGR_018', [5; -25] %[30; -10]
            'MGR_019', [20; -5]
            'MGR_028', [20; -5] %[40; 5]
            'MGR_033', [30; 0] %[55 20]
            };
end
rowid = ~cellfun(@isempty,strfind(map(:,1),name));
if all(rowid==0)
    fprintf('This study [%s] / joint combination has not been recorded\n',name)
    keyboard
end
off = map{rowid, 2};

end

% ------------------------------------------------------------------------
function T = selectRefSlices(X,plane)
T = cell(size(X));
for p = 1:numel(X)
    D = X{p};
    % Find which slice we want by working out which one has the minimum
    % distance to the action plane
    d = NaN(1,numel(D));
    for s = 1:numel(D)
        ds = distancePointPlane(D{s},plane);
        d(s) = sqrt(sum(ds.^2)) / numel(ds);
    end
    [~,s] = min(abs(d));
    T{p} = D{s};
end

end %selectRefSlice()



% ------------------------------------------------------------------------
function X = trunc(X,plane,opt)
% OPT:
%   below - discard points below plane
%   above - discard points above plane

% Cell drilling:
if iscell(X)
    for j = 1:numel(X)
        X{j} = trunc(X{j},plane,opt);
    end
else
    % Main selection function:
    below = isBelowPlane(X,plane);
    switch opt
        case 'below'
            % Keep everything that is above the plane:
            X = X(~below,:);
        case 'above'
            % Keep everything that is below the plane
            X = X(below,:);
    end
end

end %trunc()


% ------------------------------------------------------------------------
function X = all_tendon_data(model)
X = {model.LoRes.xyz}';  % That's all our data.

%
% Get trace plane normal:
n = get_slice_plane_normal(cat(1,X{:}));

% Now partition into slices:
for j = 1:numel(X)
    X{j} = split_cloud_to_slices(X{j},n);
end
%}

% Now to split them into slices:



% And now we're in the form we need...

end %all_tendon_data()


% ------------------------------------------------------------------------
function p = midSaggitalPlane(plane,cg,cs)
% plane - the nearest plane on which a trace is defined
%    cg - the centroid of the cloud
%    cs - principal component coordinate system 
pta = cg;
ptb = cg + cs(:,1)';
pt1 = projPointOnPlane(pta,plane);
pt2 = projPointOnPlane(ptb,plane);

v1 = pt2 - pt1;
v1 = v1/norm(v1);   %This now is the "up" vector of the trace plane

n = planeNormal(plane);

v2 = cross(v1,n);   % This now is the "forward" vector of the trace plane

% Get required tilt angle of plane:
ptb_dash = ptb + (pt1-pta);
tilt = anglePoints3d(pt2,pt1,ptb_dash); % [rad]

% Now we can rotate that plane about its second direction (v2) by angle tilt
R = createRotation3dLineAngle(v2,tilt);
R = R(1:3,1:3);
v2d = v2*R;     % actually v2d == v2 because it's the rotation axis
v1d = v1*R;

% Now plane can be re-defined with the rotated vectors:
plane(4:6) = v1d;
plane(7:9) = v2d;

% And set the hang point as the cg:
plane(1:3) = cg;

% And hand out the modified plane:
p = plane;


end %midSaggitalPlane()

% ------------------------------------------------------------------------
function p = nearest_scan_plane(xyz,pt)
%xyz is a N-by-1 cell array of N-by-3 traces (slices) for a single phase
% Define a plane by selecting 3 points on the appropriate slice:

for j = 1:numel(xyz)
    tj = xyz{j};
    p = planeThroughTrace(tj);
    d(j) = distancePointPlane(pt,p);
end
[dmin,j] = min(abs(d));
p = planeThroughTrace(xyz{j});

end %nearest_scan_plane()


% ------------------------------------------------------------------------
function T = split_cloud_to_slices(xyz,nml)
% Divide point cloud up into groups base on slice plane normal vector nml

o = mean(xyz);

% Project all points onto the plane normal & find the distance 
%xp = project_points_onto_line(xyz,repmat(n,r,1));
%ds = distancePointPlane(xyz,plane);
d = -sum(bsxfun(@times, nml, bsxfun(@minus, o, xyz)), 2);

scale3 = (max(xyz) - min(xyz));
scale1 = sqrt(sum(scale3.^2));
TOL = scale1/1e8;

[~,ia,ib] = unique(round(d/TOL)*TOL);

ns = numel(ia);
T = cell(ns,1);
for j = 1:ns
    T{j} = xyz( ib == j, :);
end


end 

% ------------------------------------------------------------------------
function nml = get_slice_plane_normal(xyz)
% Assume that points in xyz are ordered.  
% This won't work if they're not.

dj = 5;
np = size(xyz,1);
for j = 4:dj:np
    [Cj,~,Lj] = princomp(xyz(1:j,:),'econ');
    if Lj(3) > eps
        break
    else
        C = Cj;
        L = Lj;
    end
end
j = j - dj;
nml = cross(C(:,1),C(:,2))';
p0 = mean(xyz(1:j,:));
p = createPlane(p0,nml);
d = distancePointPlane(xyz,p);
du = unique(round(d/1e-10)*1e-10,'rows');

if numel(du) > 8
    fprintf('Problem: we shouldn''t have more than 8 slices, normally.\n');
    keyboard
end

end


% ------------------------------------------------------------------------
function n = get_slice_plane_normal_old(xyz)
% Find the normal vector of the slice plane that is used to define the
% traces that constitute the N-by-3 cloud xyz
%
% We can do this by repeatedly taking four points at random. If all
% four points are coplanar, we have a candidate plane.  Store its normal,
% repeat a million times, and find the normal that is repeated (with
% tolerance)
error('This function isn''t reliable -- doesn''t always solve correctly!!')
r  = size(xyz,1); % # rows
o  = mean(xyz);
% Basic vector projection is like this:
% project = @(a,b) dot(a,b,2)*b/norm(b);

% In matrix form, it looks like this:
% bsxfun(@rdivide, bsxfun(@times,dot(a,b,2),b), norm(b) )

% But we also have to define @norm in the conventional way:
norm = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);
tic
N = [];%NaN(ni,3);
j = 1;
nrepeats = 5;
while toc < 30  % max running time of 30s 
    % Select some indices to define a random plane:
    idx = select_indices('random',3);
    pts = xyz(idx,:);
    
    % Check that no 3 points are colinear:
    combs = nchoosek(1:numel(idx),3);
    for k = 1:size(combs,1)
        if iscolinear( pts(combs(1,:)',:) )
            disp('Skipping set with reduced dimensionality')
            continue
        end
    end
    
    inds = [1:size(xyz,1)]';
    inds(idx) = [];
    xyz_red = xyz(inds,:);
    
    
    d = distancePointPlane(xyz_red,createPlane(pts(1,:),pts(2,:),pts(3,:)));
    np = sum( abs(d) < eps );
    
    if np == 0
        continue
    end
    
    % So now we have at least 4 coplanar points, which makes for a candidate plane;
    % calculate & store normal
    
    % Two vectors in plane:
    u = pts(3,:) - pts(2,:); u = u./norm(u);
    v = pts(1,:) - pts(2,:); v = v./norm(v);
    
    % Plane normal:
    n = cross(u,v); n = n./norm(n); % vector norm
    
    % Try to avoid parallel vectors with opposite sign; set sign so that
    % largest absolute value component takes a positive sign
    [~,i] = max(abs(n));
    n = n*sign(n(i));
    
    % Store:
    N(j,:) = round(n*1e8)/1e8;
    
    % STOP CONDITION
    %    We want to see the normal vector repeated at least 5 times; then
    %    we can be satisfied that we've found it.
    if j >= nrepeats
        [C, ia, ib] = unique(N, 'rows');
        if j - numel(ia) + 1 >= nrepeats
            break
        end
    end
    j = j+1;
end

% Now our normal the most frequently repeated element
n = C(mode(ib),:);


    % ---------------------------------
    function inds = select_indices(type,n)
        switch type
            case 'random'
                % randomly select 3 unique indices:
                inds = round( rand(10,1) * (r-1)) + 1;
                inds = unique(inds);
                inds = inds(1:n);
            case 'predict1'
            case 'predict2'
        end
        
    end

end %get_slice_plane_normal()


% ------------------------------------------------------------------------
function xp = project_points_onto_line(pts,vec)
xp = bsxfun(@rdivide, bsxfun(@times,dot(pts,vec,2),vec), norm(vec) );
end


% ------------------------------------------------------------------------
function tf = iscolinear(p1,p2,p3)
if nargin == 1
    pts = p1;
    p1 = pts(1,:);
    p2 = pts(2,:);
    p3 = pts(3,:);
end
%term1 =  abs(dot(cross(p1, p2),p3)) ;
term =  norm(cross(p2-p1, p3-p1)) ;
%keyboard
tf = term < eps;
end
function theta = joint_dt_proxy(mdls)
% Get joint angle proxy


% Create a proxy for foot angle, used in smoothing:
%theta = foot_angle_proxy(mdls);
%theta = centroid_angle_proxy(mdls);
%theta = in_plane_angle_proxy(mdls);

theta = step_proxy_from_helical();
%theta = smooth(theta,0.5,'lowess');
%theta = 1:22;


% ------------------------------------------------------------------------
function theta = step_proxy_from_helical()
% This is such a mess....

handles = guidata(findall(0,'type','figure','Name','Registration'));

%% Get helical axis & rotation angle
hax = calcHelicalAxes(handles.HelicalAxis,handles.Models,'raw');
hax = hax(1); % For safety
qh = quaternion.AngleAxis(hax.Angle,hax.Axis);
[dtheta,axs] = qh.angleaxis;    % now we have dtheta with sign?
axs = axs';

% Unwrap
% a = zeros(size(dtheta));
% for j = 2:numel(dtheta)
%     a(j) = vectorAngle3d(axs(j-1,:),axs(j,:));
% end
% flip = a>(pi/2);
% axs(flip,:) = -axs(flip,:);
% dtheta(flip) = -dtheta(flip);

% flip = false(size(dtheta));
% for j = 2:numel(dtheta)
%     a = vectorAngle3d(axs(j-1,:),axs(j,:));
%     if a > pi/2
%         flip(j:end) = ~flip(j:end);
%     end
% end
% axs(flip,:) = -axs(flip,:);
% dtheta(flip) = -dtheta(flip);

%% Now use only the in-plane component
names = {handles.Models.Tag};
m1 = handles.Models(strcmpi(names,handles.HelicalAxis(1).Item1));
m2 = handles.Models(strcmpi(names,handles.HelicalAxis(1).Item2));

%{
% create a mid-sagittal plane using all the dynamic traces from the first
% phase:  (hack)
X = [];
for j = 1:numel(handles.Models)
    X = [X; handles.Models(j).LoRes(1).xyz];
end
%}
np = numel(m1.q);
% Here we find the traces of all the centroids, and then put a plane
% through the those points - closest thing to finding the plane of action.
for j = np:-1:1
    q1 = m1.q(j);
    x1 = m1.x(j,:);
    q2 = m2.q(j);
    x2 = m2.x(j,:);
    %CS1 = q1.rotate(m1.HiRes.CS);
    %CS2 = q2.rotate(m2.HiRes.CS);
    %qrel(j,1) = quaternion.RotationMatrix( CS2'*CS1 );
    % Get trace of origins, for later
    o1(j,:) = q1.rotate(m1.HiRes.Origin) + x1; % rotating (0,0,0) is unnecessary, but whatever.
    o2(j,:) = q2.rotate(m2.HiRes.Origin) + x2; %   "   "
end
X = [o1;o2];

% Now create plane using pca:
PC = princomp(X);
u = PC(:,1)';
v = PC(:,2)';
nml = cross(u,v);

% Use only in-plane component:
alfa = vectorAngle3d(axs,nml);
dtheta = dtheta.*cos(alfa);

% cumsum can't handle NaNs
dtheta(isnan(dtheta)) = 0;
theta = cumsum( (abs(dtheta)) );        % Abs? not abs?

%keyboard

% ------------------------------------------------------------------------
function theta = in_plane_angle_proxy(mdls)
% Find the angle between the first principal components 

names = {mdls.Tag}; 
sel = [];
while numel(sel) ~= 2
    [sel,ok] = listdlg('ListString',names,...
        'SelectionMode','multiple',...
        'ListSize',[180 (numel(names)+1)*13],...
        'PromptString','Select the two models for comparison', ...
        'Name','Select models' ...
        );
    if ~ok
        error('User cancelled')
    end
end

np = numel(mdls(1).q);
% Here we essentially do calculation to find the IHA angle & direction
for j = np:-1:1
    q1 = mdls(sel(1)).q(j);
    x1 = mdls(sel(1)).x(j,:);
    q2 = mdls(sel(2)).q(j);
    x2 = mdls(sel(2)).x(j,:);
    CS1 = q1.rotate(mdls(sel(1)).HiRes.CS);
    CS2 = q2.rotate(mdls(sel(2)).HiRes.CS);
    qrel(j,1) = quaternion.RotationMatrix( CS2'*CS1 );
    % Get trace of origins, for later
    o1(j,:) = q1.rotate(mdls(sel(1)).HiRes.Origin) + x1; % rotating (0,0,0) is unnecessary, but whatever.
    o2(j,:) = q2.rotate(mdls(sel(2)).HiRes.Origin) + x2; %   "   "
end
% Now theta is the angle required to rotate mdl(sel(1)) to align with
% mdl(sel(2))
[theta,axs] = qrel.unwrap.angleaxis;
axs = axs';


% Then from there we find a mean plane through the centroid traces
% Alternatively, we could put a mean plane through all the dynamic traces
% (so long as the tendon is included) and that might be good too.
PC = princomp([o1;o2]);
u = PC(:,1)';
v = PC(:,2)';
nml = cross(u,v);

% Angle between axs and plane normal:
alfa = vectorAngle3d(axs,nml); 
%keyboard
% Now factor theta by sin(alfa) to the get amount of in-plane rotation:
theta = theta.*sin(alfa);


% ------------------------------------------------------------------------
function theta = centroid_angle_proxy(mdls)
% Use positions of centroids to calculate a proxy for the angle.

names = {mdls.Tag}; 

% Find the model through which the axis runs.
search_for = @(objname) ~cellfun(@isempty,strfind(lower(names),objname));
axsId = search_for('tal');

if ~any(axsId)
    axsId = search_for('axis');
end

if sum(axsId)~=1
    [sel,ok] = listdlg('ListString',names,...
        'SelectionMode','single',...
        'ListSize',[180 (numel(names)+1)*13],...
        ...%'InitialValue',1,...
        'PromptString',{'Axis model could not be identified.','',...
        'Please select the model through which the axis runs'}, ...
        'Name','Select angle reference' ...
        );
    if ~ok
        error('User cancelled')
    end
    axs = mdls(sel);
else
    axs = mdls(axsId);
end

other = mdls(~axsId);

obj1 = other(1);
obj2 = other(2);

addpath(genpath('tendon_loa_wip'))

p1 = reshape([obj1.LoRes.Origin],3,[])';
p2 = reshape([axs.LoRes.Origin],3,[])';
p3 = reshape([obj2.LoRes.Origin],3,[])';

l1 = createLine3d(p1,p2);
l2 = createLine3d(p2,p3);

theta = vectorAngle3d(l1(:,4:6),l2(:,4:6));


% The angle must have a sign for plantarflexion / dorsiflexion, so we need
% to project it into 2D
PC = princomp([p1;p2;p3]);
o = mean([p1;p2;p3]);
rp = createPlane(o,PC(:,3)');   % Mean 2D rotation plane
pn = planeNormal(rp);
rv = normalizeVector3d(cross(normalizeVector3d(l1(:,4:6)),normalizeVector3d(l2(:,4:6))));

% Now factor theta by the component in-plane, and the sign:
scale = vectorNorm(project(rv,pn));
pos = bsxfun(@minus,rv,pn);
neg = bsxfun(@minus,rv,-pn);
sgn = (vectorNorm(pos) < vectorNorm(neg))*2-1;

theta = theta.*scale.*sgn;




%{

drawpoint = @(a,varargin)plot3(a(1),a(2),a(3),varargin{:});
drawline = @(a,b,varargin)plot3([a(1);b(1)],[a(2);b(2)],[a(3),b(3)],varargin{:});
%%
figure,
for j = 1:size(l1,1)
    drawpoint(p1(j,:),'r*')
    hold on
    drawpoint(p2(j,:),'k*')
    drawpoint(p3(j,:),'b*')
    
    drawline(p1(j,:),p2(j,:),'r');
    drawline(p2(j,:),p3(j,:),'b');
    
    pn = p2(j,:) + n(j,:)*30;
    
    
    drawline(p2(j,:),pn,'g')
    
    grid on
    axis equal
    axis auto
    %view([90,0])
    pause
    hold off
end
%}
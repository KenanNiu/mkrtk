% ------------------------------------------------------------------------
function proxy = foot_angle_proxy(mdls)
% Create a proxy for foot angle, based off the orientation of the calcaneus

names = {mdls.Tag}; 
calcId = ~cellfun(@isempty,strfind(lower(names),'calc'));
if all(calcId==0)
    [sel,ok] = listdlg('ListString',names,...
        'SelectionMode','single',...
        'ListSize',[180 (numel(names)+1)*13],...
        ...%'InitialValue',1,...
        'PromptString',{'Calcaneus could not be identified.','',...
        'Please select the model to use as a proxy for the foot angle'}, ...
        'Name','Select angle reference' ...
        );
    if ~ok
        error('User cancelled')
    end
    calc = mdls(sel);
else
    calc = mdls(calcId);
end


% Action plane bisects the entire joint into roughly medial/lateral sides:
actionplane = get_action_plane(mdls);
nml = planeNormal(actionplane);

% Now our joint angle proxy is going to be the projection of the first
% principal component of the calcaneus onto the action plane

if isfield(calc,'qraw')
    q = calc.qraw;
else
    q = calc.q;
end

vhat = [1 0 0];
n = numel(q);
proxy = NaN(n,1);
for j = 1:n
    vj = q(j).rotate(vhat);                 % Find orientation of first principal component vector
    vj = project_vector_on_plane(vj,actionplane);   % Project vector onto plane
    if j == 1;
        v1 = vj;
    end
    proxy(j) = signed_vectorAngle3d(nml,v1,vj);
    
end


% ------------------------------------------------------------------------
function actionplane = get_action_plane(mdls)
% Let's take all the models in position 1, and find the first two principal
% components - let's call this our median plane
phase = 1;
n = numel(mdls);
C = cell(n,1);
for j = 1:n
    mj = mdls(j);
    if isempty(mj.HiRes)
        xyz = [];
    else
        cld = mj.HiRes.transform(mj.q(phase),mj.x(phase,:));
        xyz = cld.xyz;
    end
    C{j} = xyz;
end
XYZ = cat(1,C{:});

% Now get principal components:
coeff = princomp(XYZ);
coeff(:,3) = cross(coeff(:,1),coeff(:,2)); % This ensures a righ-hand axis system, since princomp doesn't care

% And use that to define a plane:
cg = mean(XYZ);
actionplane = [cg coeff(:,1)' coeff(:,2)'];



% ------------------------------------------------------------------------
function v2 = project_vector_on_plane(v,p)
% use geom3d's projPointOnPlane function

pt0 = [0 0 0];
pt1 = v;

pt0d = projPointOnPlane(pt0,p);
pt1d = projPointOnPlane(pt1,p);

v2 = pt1d-pt0d;

v2 = v2/ sqrt( sum(v2.^2) );    % normalize


% ------------------------------------------------------------------------
function theta = signed_vectorAngle3d(nml,v1,v2)
% Calculate the signed angle that v1 must be rotated around nml to arrive
% at v2

tol = 1e-5;
isEqualWithTol = @(a,b) isequal(size(a),size(b)) && all( abs( a-b ) < tol);
normalize = @(v)v./sqrt( sum(v.^2) );

nml = normalize(nml);
v1  = normalize(v1);
v2  = normalize(v2);

cp = normalize(cross(v1,v2));

theta = vectorAngle3d(v1,v2);
%subspace(v1',v2') % same thing

if theta == 0  % vectors are parallel, cross(v1,v2)==[0 0 0]
    sgn = +1;
elseif isEqualWithTol( nml, cp )
    sgn = +1;
elseif isEqualWithTol( nml, -cp )
    sgn = -1;
else
    warning('nml is not a mutual normal')
    keyboard
end

theta = theta*sgn;


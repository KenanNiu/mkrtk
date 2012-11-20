function bone_class_spec

% Class specification & tests for Bone class.
%
% This class supersedes the "Models" structure.


% ======================== Basic Listing ================================ %


% Constructors
b = Bone();


% Properties
b.Tag;
b.HiRes;
b.LoRes;
b.smoothed;

% Dependent properties
b.q;     % gets qsmooth if it's not empty, otherwise gets qraw
b.x;     % same for xsmooth

% Hidden properties
b.qraw;
b.xraw;
b.qsmooth;
b.xsmooth;

% Hidden & protected properties:
b.Version;


% Methods (logical returns)
%b.isregistered;


% Methods (value returns)


% Methods (object manipulators)
b.smoothpose(@(y)y);
%b.clearpose; or clearsmoothing? or resetsmoothing


%% ======================== Variables ==================================== %

[c1s,c2s,c3s, c1d,c2d,c3d, q1,q2,q3, x1,x2,x3] = create_vars();

tags = {'First','Second','Third'};

%% ======================== Formal Testing =============================== %



% -------------- Constructors ----------------- %
bstruct = struct(...
    'Tag',tags,...
    'HiRes',{c1s,c2s,c3s},...
    'LoRes',{c1d,c2d,c3d},...
    'q',{q1,q2,q3},...
    'x',{x1,x2,x3},...
    'qraw',{q1,q2,q3},...
    'xraw',{x1,x2,x3} );
bstruct = bstruct(:);

b = Bone(bstruct);      % Create from struct (backward compatibility)

b = Bone();             % Blank bone object
b = repmat(b,3,1);      %  ...assign properties below
%%

% -------------- Assign Properties ----------------- %
[b.HiRes] = deal(c1s,c2s,c3s);
[b.LoRes] = deal(c1d,c2d,c3d);

[b.q] = deal(q1,q2,q3);
[b.x] = deal(x1,x2,x3);
[b.Tag] = deal(tags{:});


% -------------- Method Calls ----------------- %
theta = [1:numel(q1)];
%span  = 0.5;
%fun = @(y)smooth(theta,y,span,'rloess');
span = numel(q1)-1;
fun = @(y)smooth(theta,y,span,'moving');

assertFalse(any([b.smoothed]))

b = b.smoothpose(fun);

assertTrue(all([b.smoothed]))

b = b.clearsmoothing;

assertFalse(any([b.smoothed]))




% -------------- Display ----------------- %
hf = figure;
ax = axes;
grid(ax,'on')
hold(ax,'on')
clrs = lines(numel(b));
for p = 1:numel(q1)
    for bj = 1:numel(b)
        b(bj).LoRes(p).plot(ax,'color','k','Marker','*')
        b(bj).HiRes.transform(b(bj).q(p),b(bj).x(p,:)).plot(ax,'Color',clrs(bj,:),'Marker','.','MarkerSize',7)
    end
    axis(ax,'tight','equal')
    if p==1, view(3), end
    pause
    cla(ax)
end
close(hf);


end



%=========================================================================
function [c1s,c2s,c3s, c1d,c2d,c3d, q1,q2,q3, x1,x2,x3] = create_vars() %#ok<STOUT>

NB = 3;     % Number of bones
ND = 10;    % Number of dynamic positions
% Create states with which to perturb the dynamic clouds:
%  - q1, x1
%  - q2, x2
%  - q3, x3
q1 = quaternion.AngleAxis(0 + linspace(0,pi/8,ND)', repmat([1 0 0],ND,1));      
q2 = quaternion.AngleAxis(pi/4 + linspace(0,pi/8,ND)', repmat([0 1 0],ND,1));   
q3 = quaternion.AngleAxis(repmat(-pi/8,ND,1), normr(interp1([0 1 0; 1 0 0],linspace(1,2,ND))));

x1 = interp1( [0 0 0; 0 0 0], linspace(1,2,ND) );       
x2 = interp1( [0 90 90; 0 120 120], linspace(1,2,ND) );   
x3 = interp1( [-40 0 0; -80 0 0], linspace(1,2,ND) );   

% Scale:
S = [100, 50, 80];

% Create Clouds:
for j = 1:NB
    % Create static clouds:
    %  - c1s
    %  - c2s
    %  - c3s
    Z = membrane(j,90);   % Static
    [X,Y] = meshgrid(1:size(Z,1),1:size(Z,2));
    f = max([max(X(:)) max(Y(:))]); %\
    X = X./f*S(j);                  % \_ Scale
    Y = Y./f*S(j);                  % /
    Z = Z*S(j);                     %/
    XYZ_s = [X(:) Y(:) Z(:)]; %#ok<NASGU>
    % eval produces (j==3):
    %   c3s = Cloud(XYZ_s);
    eval(sprintf('c%ds = Cloud(XYZ_s);',j));
    
    % Create dynamic clouds:
    %  - c1d
    %  - c2d
    %  - c3d
    Z = membrane(j,10); % Dynamic
    [X,Y] = meshgrid(1:size(Z,1),1:size(Z,2));
    f = max([max(X(:)) max(Y(:))]); %\
    X = X./f*S(j);                  % \_ Scale
    Y = Y./f*S(j);                  % /
    Z = Z*S(j);                     %/
    XYZ_dk = [X(:) Y(:) Z(:)];
    cd0 = Cloud(XYZ_dk); %#ok<NASGU>
    
    for k = ND:-1:1
        % eval produces (j==3):
        %   c3d(k,1) = cd0.transform(q3(k),x3(k,:));
        eval(sprintf('c%dd(k,1) = cd0.transform(q%d(k),x%d(k,:));',j,j,j));
    end
end

end %create_vars()


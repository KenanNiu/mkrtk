function cloud_class_spec
% CLOUD Class spec & tests for CLOUD class



%#ok<*NASGU>

% Preparation
assertEqual(exist('Cloud.m','file'),2)
assertEqual(isa(Cloud(),'Cloud'),true)

close all

% ========================== Variables ================================== %

% Basic data:
xyzNx3 = rand(50,3);    
xyzNx3xN = rand(50,3,3);

% Formatted as traces:
X = {rand(10,3); 
    rand(50,3); 
    rand(43,3); 
    rand(60,3)};

% Source data files:
here = fileparts(mfilename('fullpath'));
dfname = 'MGR_033_ankle_talus_dynamic.mat';
sfname = 'MGR_033_ankle_tibia_static.mat';
dynamicfile = [here filesep 'cloud_data' filesep dfname];
staticfile  = [here filesep 'cloud_data' filesep sfname];

% Scalar transformation variables:
t = [10,-5,30];             % Translation vector
ang = pi/2;                 %\_ Angle-axis rotation
axs = [1 0 0];              %/
q = quaternion.AngleAxis(ang,axs);  % equivalent quaternion
R = q.rotationmatrix;               % equivalent rotation matrix
rctr = [3,4,5];             % Rotation centre

% Vectorised transformation variables:
nV   = 7;
tV   = [ (linspace(0,10,nV))', (linspace(5,0,nV))', (linspace(20,30,nV))' ];
angV = (linspace(0,pi/2,nV))';  % Vary the angle
axsV = repmat(axs,nV,1);        % Copy
qV   = (quaternion.AngleAxis(angV,axsV));
RV   = qV(:).rotationmatrix;       % Convert->3x3xN
rctrV= repmat(rctr,nV,1);       % Copy

% ======================== Basic Listing ================================ %

% Constructors
c = Cloud(xyzNx3);           % [N-by-3]
c = Cloud(xyzNx3');          % [3-by-N]
c = Cloud(X);               % {M-by-1 cell array}
c = Cloud(X');              % {1-by-M cell array}

c = Cloud.Load(staticfile); % Load from file

c = c(ones(nV,1));          % To test methods for vectorisation

% Configuration:
c.setorigin;                % Set/reset default origin to mean(c.xyz)
c.setprincompcs;            % Set/reset coordinate system to default principal component CS

% These only work on single objects:
c(1).traces;                                % Return traces (if loaded from a set of ROIs or cell array)
c(1).downsample(int32(1000));               % Downsample to specified number of points.
c(1).downsample(false(size(c(1).xyz,1),1)); % Downsample with binary selector

% Translations
c.translate(t);                 % Translate multiple clouds with same vector
c.translate(tV);                % Translate multiple clouds with different vectors
c(1).translate(tV);             % Translate one cloud into multiple locations


% Rotations - Rotation matrix
c(1).rotate(R);                 % Rotate single cloud with matrix
c(1).rotate(R,rctr)             %   '-> around origin

c(1).rotate(RV);                % Rotate single cloud into multiple positions
c(1).rotate(RV,rctrV);          %   '-> around origin

c.rotate(R);                    % Rotate multiple clouds with same matrix
c.rotate(R,rctrV);              %   '-> around origin

c.rotate(RV);                   % Rotate multiple clouds with separate matrices
c.rotate(RV,rctrV);             %   '-> around origin


% Rotations - Angle/Axis
c(1).rotate(ang,axs);           % Rotate single cloud with single angle/axis
c(1).rotate(ang,axs,rctr);      %   '-> around origin

c(1).rotate(angV,axsV);         % Rotate single cloud into multiple positions
c(1).rotate(angV,axsV,rctrV);   %   '-> around origin

c.rotate(ang,axs);              % Rotate multiple clouds with same angle/axis
c.rotate(ang,axs,rctr);         %   '-> around origin

c.rotate(angV,axsV);            % Rotate multiple clouds with separate angle/axis
c.rotate(angV,axsV,rctrV);      %   '-> around origin


% Rotations - Quaternion
c(1).rotate(q);                 % Rotate single cloud with quaternion
c(1).rotate(q,rctr)             %   '-> around origin

c(1).rotate(qV);                % Rotate single cloud into multiple positions
c(1).rotate(qV,rctrV);          %   '-> around origin

c.rotate(q);                    % Rotate multiple clouds with same quaternion
c.rotate(q,rctrV);              %   '-> around origin

c.rotate(qV);                   % Rotate multiple clouds with separate quaternions
c.rotate(qV,rctrV);             %   '-> around origin


% Transformations: Rotation then translation
c(1).transform(R,t);            % Transform single cloud with single tranformation
c(1).transform(RV,tV);          % Transform single cloud into mulitple positions

c.transform(R,t);               % Transform multiple clouds with single transformation
c.transform(RV,tV);             % Transform multiple clouds with separate transformations



% Recovering transformations

gettransform(c,c,'mat');        % Return rotation matrices & translation vector: [R,t]
gettransform(c,c,'quat');       % Return quaternions & translation vectors: [q,t]


% Plotting / graphics
c.plot;                         % Plot cloud
c.plot(gca);                    % Plot to specified axes
c.plot('r', 'Marker','*');      % With options
c.plotcs;                       % Plot coordinate system


clear c


% ======================== Formal Testing =============================== %
close all



% -------------- Constructors ----------------- %

% Loading from file
cs = Cloud.Load(staticfile );
cd = Cloud.Load(dynamicfile);

assertEqual( numel(cs), 1 )
assertEqual( numel(cd), 33 )

% Normal constructors:
assertEqual( 1, numel( Cloud() ) )      % 1-by-1 with empty fields
assertEqual( 0, numel( Cloud([]) ) )    % 0-by-0
assertEqual( 0, numel( Cloud.empty ) )  % 0-by-0 (prefereable to Cloud([]) )

% Check cloud creation from matrices:
% N-by-3
c = Cloud(xyzNx3);
assertEqual(xyzNx3, c.xyz)
% 3-by-N
c = Cloud(xyzNx3');
assertEqual(xyzNx3, c.xyz)

% Failing cases
assertExceptionThrown( @()Cloud(xyzNx3xN), '' )     % too many dimensions
assertExceptionThrown( @()Cloud(rand(20,2)), '' )   % too few columns


% Cloud creation from cell:
cX = Cloud(X);
cx = Cloud(cat(1,X{:}));

assertEqual(cX.xyz,cx.xyz)

% Converting a structure (version upgrading)
    s.Tag = 'Test';
    s.xyz = xyzNx3;
    s.Origin = [5.1 9.8 -7.2];
cS = Cloud(s);
    
    wstate = warning('off','MATLAB:structOnObject');
    s = (arrayfun(@struct,cd))';
    
cS = Cloud(s);

    warning(wstate)

% Just check a few properties:
for j = 1:numel(cS)
    assertEqual( cS(j).CS,  s(j).CS  );
    assertEqual( cS(j).Tag, s(j).Tag );
    assertEqual( cS(j).xyz, s(j).xyz );
    assertEqual( cS(j).Origin, s(j).Origin );
end



% -------------- Instance Methods ----------------- %

% Traces:
assertEqual(size(cs.traces),    [53 1])     % cs.traces
assertEqual(size(cd(1).traces), [ 7 1])     % cd(1).traces


% -------------- Plotting ----------------- %
    hf = figure;
    
% Simple plot
hc = cs.plot;        

    clf(hf)
    ax = axes;
    hold(ax,'on')
    grid(ax,'on')

% Specify axes
cs.plot(ax);        
cd(1).plot(ax,'LineWidth',2,'Marker','+','color','k'); % Specify properties


    axis(ax,'tight','equal')
    view(3)
    
% Plot coordinate system  
hcs = cs.plotcs;    
    
    hold(ax,'off')
    
% Plot traces
cd.plottraces(ax)


% -------------- Performing Transformations ----------------- %
fprinf(2,'This needs to be updated to check vectorised solutions...\n');
c = cs;
% Translation:
ct = c.translate(t); 
assertAlmostEqual(mean(c.xyz)+t, mean(ct.xyz))

    figure, c.plot('b'), grid on, axis equal
    
% Rotation
cr = c.rotate(ang,axs,c.Origin);
cr = c.rotate(R);

    hold on
    
cr.plot('color',[0 0.7 0])
cr.plotcs;
    
    axis tight

% Rotation + Translation
cx = c.transform(q,t);

% -------------- Recovering Transformations ----------------- %
[Rx,tx] = gettransform(c,cx,'mat');
[qx,tx] = gettransform(c,cx,'quat');
assertAlmostEqual(R, Rx, 1e-12)
assertAlmostEqual(q, qx, 1e-12)
assertAlmostEqual(t, tx, 1e-12)


% -------------- Down-sampling & Selection ----------------- %
% Down-sampling to N points:
    np = int32(1000);
    
cs_coarsen  = cs.downsample(np);   % integer => # of points
%cs_coarseds = cs.downsample( 1.5 );         % double  => min spacing

dsplot = @(cld)plot(cld,'r','Marker','*','LineStyle','none');
figure, 
plot(cs,'b-');
hold on, grid on, axis equal, axis tight
hc = dsplot(cs); 
title('Recursive downsampling test')
while size(cs.xyz,1) > 10
    pause(0.7)
    np = size(cs.xyz,1);
    cs = cs.downsample( int32(np/2) );
    delete(hc);
    hc = dsplot(cs);
end
close all

% Down-sampling by selection
    keep = false(size(cs.xyz,1),1);
    keep(1:3:size(cs.xyz)) = true;
cs_coarsei = cs.downsample(keep);

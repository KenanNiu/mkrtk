clear all

% Geometry classes:
%   - Point:   (x,y,z)
%   - Vector: (xhat, yhat, zhat)
%   - Line :   ( point, vector )  OR segment: ( point, point )
%   - Plane:   ( point, normal )  OR: ( point, xhat, yhat )
%

% Can they all inherit properties of "GEOMETRY" superlass? (maybe GEOMETRY3D)

% LINE3D thoughts:
%   - Two basic forms:
%       - ray
%       - segment
%   - Can we superclass? line3d would then be abstract:
%       - ray << line3d
%       - segment << line3d
%       - polyline << segment
%   - If we superclass, need to work out which methods go in the
%     superclass...
%
% Operations
%   - Vector -> normalize
%            -> define cross product & dot product
%   - intersect -> with line or plane
%   - project
%   - orthogonal distance -> with another line -> produces a point-point line as the result 
%   - angle between

% Angle between two vectors ->  SUBSPACE()
%  (also angle between vector & plane?  or dihedral angle of two planes?)



initTestSuite;

% ------------------------------------------------------------------------
% NOTES
%   - Actions that act on two different classes should be handled by the
%     higher class.  Ie, 
%           - projecting a point on a line  -> line class
%           - projecting line on plane      -> plane class


%superclass?? --> geom3d


% Methods:
%   - project()
%   - rotate()
%

% ------------------------------------------------------------------------
% POINT

translation = [-10 4 -5];
scale = 10;
a_p = [0 3 4];
b_p = [7 8 9];
c_p = a_p + translation;
d = a_p * scale;

[dim1, dim2, dim3] = deal(10, 6, 8);

A  = repmat(a_p', [1, dim1, dim2, dim3]);
A1 = repmat(a_p,  [dim1, 1, dim2, dim3]);
B  = repmat(b_p', [1, dim1, dim2, dim3]);
C  = repmat(c_p', [1, dim1, dim2, dim3]);
D  = repmat(d', [1, dim1, dim2, dim3]);

% Basic constructors:
assertTrue( isa(point3d(), 'point3d') );
pa = point3d(a_p);
pb = point3d(b_p);
pc = point3d(c_p);
pd = point3d(d);

assertEqual( point3d(a_p'), pa )

% Matrix constructors:
pA1 = point3d(A1);
pA = point3d(A);
pB = point3d(B);
pC = point3d(C);
pD = point3d(D);

assertEqual( size(pA), [dim1 dim2 dim3] )
assertEqual( pA(1),    pa )
assertEqual( pA(end),  pa )

assertEqual( pA, pA1)


% Deconstructors
assertEqual( pA.double, A )
assertEqual( pA.x, squeeze(A(1,:,:,:)) )
assertEqual( pA.y, squeeze(A(2,:,:,:)) )
assertEqual( pA.z, squeeze(A(3,:,:,:)) )


% Basic equlity methods:
assertTrue(  pa == pa ) % tests eq()
assertFalse( pa == pb ) % tests eq()
assertTrue(  pa ~= pb ) % tests neq()
assertFalse( pa ~= pa ) % tests neq()

% Translation:
assertEqual( pA + translation, pC )
assertEqual( pA.translate(translation), pC )

% Scaling:
assertEqual( pA.* scale, pD )
assertEqual( pA * scale, pD )

% Angle between points
theta = angle(pa,pb,pc);
addpath(genpath('geom3d'))
if exist('anglePoints3d','file') == 2
    assertEqual( angle(pa,pb,pc), anglePoints3d(a_p,b_p,c_p) )
else
    warning('Could not run check - FEX geom3d tools not found') %#ok<WNTAG>
end
Theta_aBC = angle(pa,pB,pC);
Theta_AbC = angle(pA,pb,pC);
Theta_ABc = angle(pA,pB,pc);
assertTrue( all( Theta_aBC(:) == theta ) )
assertTrue( all( Theta_AbC(:) == theta ) )
assertTrue( all( Theta_ABc(:) == theta ) )




% ------------------------------------------------------------------------
% VECTOR

clear all

% Set up variables
norm = @(v)sqrt( sum( v.^2 ) );
normalize = @(v)v./norm(v);
a_p = [3 4 5];
b_p = [-1 -2 -3];

acrossb = cross(a_p,b_p);
adotb = dot(a_p,b_p);
theta = acos( dot( normalize(a_p), normalize(b_p) ) );

[dim1, dim2, dim3] = deal(10, 6, 8);
A  = repmat(a_p', [1,dim1,dim2,dim3]);
A1 = repmat(a_p,  [dim1,1,dim2,dim3]);
B  = repmat(b_p', [1,dim1,dim2,dim3]);
AcrossB = cross(A,B);
Theta = repmat(theta,[dim1,dim2,dim3]);


% Basic constructors:
assertTrue( isa(vector3d(),'vector3d') )
va = vector3d(a_p);
vb = vector3d(b_p);

assertEqual( va, vector3d(a_p') )

% Matrix constructors:
vA  = vector3d(A);
vA1 = vector3d(A1);
vB  = vector3d(B);

assertEqual( size(vA), [dim1, dim2, dim3] )
assertEqual( vA(1),    va )
assertEqual( vA(end),  va )

assertEqual( vA, vA1 )

% Deconstructors
assertEqual( vA.double, A )
assertEqual( vA.u, squeeze(A(1,:,:,:)) )
assertEqual( vA.v, squeeze(A(2,:,:,:)) )
assertEqual( vA.w, squeeze(A(3,:,:,:)) )


% Basic equlity methods:
assertTrue(  va == va ) % tests eq()
assertFalse( va == vb ) % tests eq()
assertTrue(  va ~= vb ) % tests neq()
assertFalse( va ~= va ) % tests neq()

% Vector scaling:
isequivalent = vA.equiv(vA*1.5);        % Test mtimes()
assertTrue( all(isequivalent(:)) )      % Test equiv() 
assertEqual( vA.*3, vector3d(A.*3) )    % Test times() (element-by-element multiplication)

% Vector addition:
assertEqual( vA + vB, vector3d(A+B) )
assertEqual( va + vB, vector3d(A+B) )
assertEqual( vA + vb, vector3d(A+B) )

% Geometric methods:
vA_norm = vA.norm;
assertEqual( vA_norm(1), norm(a_p) )
assertEqual( vA_norm(end), norm(a_p) )

vA_normalized = vA.normalize;
assertEqual( vA_normalized(1), vector3d( normalize(a_p) ) )

%   Dot & Cross products, and angles
assertEqual( cross(vA,vB), vector3d(cross(A,B)) )
assertEqual( cross(va,vB), vector3d(cross(A,B)) )
assertEqual( cross(vA,vb), vector3d(cross(A,B)) )

assertEqual( dot(vA,vB), squeeze(dot(A,B)) )
assertEqual( dot(va,vB), squeeze(dot(A,B)) )
assertEqual( dot(vA,vb), squeeze(dot(A,B)) )

assertAlmostEqual( angle(vA,vB), Theta )
assertAlmostEqual( angle(va,vB), Theta )
assertAlmostEqual( angle(vA,vb), Theta )



% ------------------------------------------------------------------------
% LINE - getting a little more complex

clear all

translation = [-10 4 -5];
scale = 10;
a_p = [0 3 4];
b_p = [7 8 9];
c_p = a_p + translation;

a_v = [3 4 5];
b_v = [1 -0.5 -.7];

[dim1, dim2, dim3] = deal(10, 6, 8);

A  = repmat(a_p', [1, dim1, dim2, dim3]);   % Point matrices
B  = repmat(b_p', [1, dim1, dim2, dim3]);   % 
C  = repmat(c_p', [1, dim1, dim2, dim3]);   %
A1 = repmat(a_p,  [dim1, 1, dim2, dim3]);   %<- alternate orientation

E  = repmat(a_v', [1, dim1, dim2, dim3]);   % Vector matrices
F  = repmat(b_v', [1, dim1, dim2, dim3]);   % 
E1 = repmat(a_v,  [dim1, 1, dim2, dim3]);   %<- alternation orientation

% Build components:
pa = point3d(a_p);  %
pb = point3d(b_p);  % Single
pc = point3d(c_p);  %
pA = point3d(A);  %
pB = point3d(B);  % Matrix
pC = point3d(C);  %

va = vector3d(a_v);   % Single
vb = vector3d(b_v);   % 
vA = vector3d(E);  % Matrix
vB = vector3d(F);  %


% ----- Line3d Class Tests:

% Existence:
assertTrue( isa(line3d(),'line3d') )
%
% (1.1) Constructor using 2 objects: point3d & vector3d
la = line3d(pa,va);
assertEqual( la.point, pa )
assertEqual( la.vector, va.normalize )

% (1.2) Constructor using 2 objects: 2x point3d
lb = line3d(pa,pb);
assertEqual( lb.point1, pa )
assertEqual( lb.point2, pb )

% (1.3) Constructor using 2 sets of doubles: point & vector
lad = line3d( 'point-vector', a_p, a_v );
assertEqual( lad.point, pa )
assertEqual( lad.vector, va.normalize )

% (1.4) Constructor using 2 sets of doubles: 2x point
lbd = line3d( 'point-point', a_p, b_p );
assertEqual( lbd.point1, pa )
assertEqual( lbd.point2, pb )
%}


% (2.1) Matrix object constructor: point3d & vector3d
lA = line3d( pA, vA );
assertEqual( size(lA), size(pA) )
assertEqual( lA.point,  pA )
assertEqual( lA.vector, vA.normalize )

% (2.2) Matrix object constructor: 2x point3d
lB = line3d( pA, pB );
assertEqual( size(lB), size(pA) )

% (2.3) Matrix double constructor: point & vector
lAd = line3d( 'point-vector', A, E );
assertEqual( size(lA), size(lAd) )

% (2.4) Matrix double constructor: 2x point
lBd = line3d( 'point-point', A, B );
assertEqual( size(lB), size(lBd) )


% Create a mixed-mode array for use later:
inds = unique(randi([1 numel(lA)],round(numel(lA)/3),1));
lM = lA; lM(inds) = lB(inds);


% Check finite-ness:
ray_A = lA.isray;
ray_B = lB.isray;
seg_A = lA.issegment;
seg_B = lB.issegment;
assertTrue( all( ray_A(:)) )
assertTrue( all(~ray_B(:)) )
assertTrue( all(~seg_A(:)) )
assertTrue( all( seg_B(:)) )
%assertTrue( ~la.isfinite )
%assertTrue(  lb.isfinite )

%assertTrue(  la.toPointVector.isfinite )


% Equality testing:
lAeq  = lA == lA;           % Test "==" works
lAdeq = lA == lAd;          % Check creation by double was the same
laeq1 = lA == la;           % Check against scalar
laeq2 = la == lA;           % ...scalar first
lAmeq = lA.eq(lA);          % Call object method directly
lBeq  = lB == lB;       % Check with p-p mode
lbeq1 = lB == lb;       % Against scalar
lbeq2 = lb == lB;       % ...scalar first
lABeq = lA == lB;       % This should be all false
assertTrue( all(lAeq(:))  )   % Evaluate tests
assertTrue( all(lAdeq(:)) )   % 
assertTrue( all(laeq1(:)) )   % 
assertTrue( all(laeq2(:)) )   % 
assertTrue( all(lAmeq(:)) )   % 
assertTrue( all(lBeq(:))  )     
assertTrue( all(lbeq1(:)) )
assertTrue( all(lbeq2(:)) )
assertFalse(all(lABeq(:)) )

% Inequality (it's just ~eq, so only need to test that it exists)
lABneq = lA ~= lB;
assertTrue( all(lABneq(:)) )


% Test getter methods:
%assertTrue( la == lad )
%assertTrue( la.eq(lad) )
assertEqual( lb, lbd )
assertEqual( la.point, lad.point )
assertEqual( la.vector, lad.vector )
assertEqual( lb.point1, lbd.point1 )
assertEqual( lb.point2, lbd.point2 )


lM.toray;
lA.toray;
lB.toray;

lM.origin;  % Just test it runs
lA.origin;
lB.origin;

lM.direction;
lA.direction;
lB.direction;


assertEqual( la.origin, la.point )
assertEqual( lb.origin, lb.point1 )




% ------------------------------------------------------------------------
% PLANE
clear all

norm = @(v)sqrt( sum( v.^2 ) );
normalize = @(v)v./norm(v);

t1 = [5  10 3];
t2 = [5 -10 3];

pt = [5 -3 8];
u  = [2 1 4];
v  = [-1 -5 -9];

uhat = normalize(u);
vhat = normalize(v);
n = cross(uhat,vhat);


[dim1, dim2, dim3] = deal(2, 5, 3);

Ptd = repmat(pt', [1, dim1, dim2, dim3] );
Ud  = repmat(u',  [1, dim1, dim2, dim3] );
Vd  = repmat(v',  [1, dim1, dim2, dim3] );
Nd  = repmat(n',  [1, dim1, dim2, dim3] );

Pt = point3d(Ptd);
U  = vector3d(Ud);
V  = vector3d(Vd);
N  = vector3d(Nd);

% Existence:
assertTrue( isa(plane3d(),'plane3d') )

% Let's not create from doubles at this stage.  It just complicates things


% ------- POINT-U-V Constructor methods: ---------- %

% From objects: point3d & vector3d & vector3d
Ppuv = plane3d( Pt, U, V );
assertEqual( size(Ppuv), size(Pt) )


% ------- POINT-NORMAL Constructor methods: ------- %

% From objects: point3d & vector3d
Ppn = plane3d( Pt, N );
assertEqual( size(Ppn), size(N) )

% ------- LINE Constructor methods: --------------- %
disp('Include these at some point')
%Plray = plane3d( line3d(Pt,N) );                    % Point-normal
%Plseg = plane3d( line3d(Pt, point3d(Ptd+Nd) ) );    % Point-point
%assertTrue( Plray == Ppn ) 
%assertTrue( Plseg == Ppn )

% ------- POINT-POINT-POINT Constructor method ------- %

% From objects: 3x point3d
Pppp = plane3d( Pt, Pt+t1, Pt+t2 );

% Option to pack points:
ptset = [point3d(pt), point3d(pt)+t1, point3d(pt)+t2];
Ppack = plane3d( ptset );


% ------- Plotting & dependencies -------- %
la = line3d('point-vector',[0 0 0],[1 1 1]);
lb = line3d('point-vector',[1 1 1],[0 1 0]);
lc = line3d('point-vector',[0 4 1],[1 2 3]);
lset = [la,lb,lc,la,lb,lc];

Ppn(1).intersect(lset);
Ppn.intersect(la);

ax = axes;
axis(ax,[-1 1 -1 1 -1 1]*50);
Ppuv(1).plot(gca,'Marker','o','MarkerFaceColor','flat','FaceColor','none')
Ppn(1).plot('m')


fprintf(2,'These are not the same -- why?\n');
% Equality testing:
eq1 = Ppuv == Ppn;
assertTrue( all(eq1(:)) )

% Check normalization:
assertTrue( plane3d(), plane3d() )




plane3d.isnormalized



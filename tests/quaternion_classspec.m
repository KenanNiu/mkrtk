clear all
initTestSuite;

% Last ported from Mark Tinknell's version of May 18th, 2012


norm = @(v)sqrt( sum( v.^2 ) );
normalize = @(v)v./norm(v);


% Some constants:
u = normalize( [1 0 0] );
v = normalize( [0 1 0] );
ax = cross(u,v);            % Mutual normal, or rotation axis
theta = acos(dot(u,v));     % Angle between u & v

x = ax(1); y = ax(2); z = ax(3);            % Calculate rotation matrix
c = cos(theta); s = sin(theta); t = (1-c);  % from axis/angle representation
R = [...                                    % see:
    t*x*x+c,   t*x*y-z*s, t*x*z+y*s         %  http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToMatrix/index.htm
    t*x*y+z*s, t*y*y+c,   t*y*z-x*s         %
    t*x*z-y*s, t*y*z+x*s, t*z*z+c   ];      %



% ------------------------------------------
% Scalar constructors

q = quaternion;

q0  = quaternion( 0, 0, 0, 0 );
qa1 = quaternion( 1, 2, 3, 4 );
qa2 = quaternion( [1 2 3 4] );
qa3 = quaternion( [1 2 3 4]');

qb = quaternion( 5 );           % [w] only, no i,j,k
qc = quaternion( 1, 2 );        % [w, i] only, no k,j
qd = quaternion( 1, 2, 3 );     % [i, j, k] only, no w


qr1 = quaternion.rand;
qr2 = quaternion.rand;

qzero = quaternion.zeros;
qnan  = quaternion.nan;
qNaN  = quaternion.NaN;

quv = quaternion.RotateUtoV( u, v );
qR  = quaternion.RotationMatrix( R );
qaa = quaternion.AngleAxis( theta, ax );

assertTrue(   q == q0  )    % initialized to zeros
assertTrue( qa1 == qa2 )    %\_ handle different input orientations
assertTrue( qa1 == qa3 )    %/

assertEqual( qb.e, [5,0,0,0]' )
assertEqual( qc.e, [1,2,0,0]' )
assertEqual( qd.e, [0,1,2,3]' )

assertFalse( qr1 == qr2    )   % must be random
assertAlmostEqual( 1, qr1.norm )   % must be normalized

assertEqual( qnan, qNaN )      % unit test uses isequalwithnans internally

assertAlmostEqual( quv.e, qR.e  ) % Check constructors
assertAlmostEqual( quv.e, qaa.e )




% ------------------------------------------
% Array constructors



% ------------------------------------------
% Conversion method tests:

[theta_aa,axs_aa] = qaa.angleaxis;
[theta_uv,axs_uv] = quv.angleaxis;
[theta_R, axs_R]  =  qR.angleaxis;

assertEqual( theta, theta_uv )
assertEqual( theta, theta_R  )
assertEqual( theta, theta_aa )


Raa = qaa.rotationmatrix;
Ruv = quv.rotationmatrix;
RR  =  qR.rotationmatrix;

assertAlmostEqual( R , Raa , 1e-15 ) 
assertAlmostEqual( R , Ruv , 1e-15 )
assertAlmostEqual( R , RR  , 1e-15 )



disp('Quaternion tests passed.')
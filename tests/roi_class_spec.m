%-------------------------------------------------------------------------
function roi_class_spec %#ok<*DEFNU>
%ROI_CLASS_SPEC Specification & tests for the 'roi' class
%
% This funtion is for the specification and testing of the ROI class using
% test-driven development through the MATLAB xUnit testing framework.
%
% See also ROI

% Basic specification - Version 2:
% --------------------------------
% r = 
%                           x: [Nx1 double]   x-vector of ROI curve
%                           y: [Nx1 double]   y-vector of ROI curve
%                       Slice: [1x1 double]   Parent slice number
%                       Phase: [1x1 double]   Parent phase number
%                       Image: [char]         Full path to image parent
%                       Color: [char]         Colorspec of the form 'r', 'red', [1 0 0], etc
%                   LineStyle: [char]         Linspec of the form '-', ':', etc
%                PixelSpacing: [1x2 double]   From DICOM header
%        ImagePositionPatient: [1x3 double]   From DICOM header
%     ImageOrientationPatient: [1x6 double]   From DICOM header
%                        Area: [-dynamic-]    Dynamically generated property
%                      Length: [-dynamic-]    Dynamically generated property
%                     Version: [1x1 double]   ROI class version number (Fixed by class)
%                         Tag: [char]         Internally generated tag, like: "ROI-42"
%
%
% Categories of tests we need to run are:
%   - Normal Constructor (from normal input arguments)
%   - Structured constructor (from struct)
%   - Version upgrading
%   - Saving / Loading
%   - Tag conficts
%   - Curve definitions (open/closed)
%   - Calculation of dependent properties
%   - Curve Adjustment (open/close, join)
%   - Conversion to 3D


%% Check we have the class on the path:
assertEqual(exist('roi.m','file'),2)    % File exists
r = roi();
assertTrue(isa(r,'roi'))                % Constructs a 'roi' object


%% Test all normal constructor methods
x = 1:10;
y = 1:10;
s = 1;
p = 1;
img = [toolboxdir('images') filesep 'imdemos' filesep 'CT-MONO2-16-ankle.dcm'];
ps = [0.5 0.5];
ipp = rand(1,3);
iop = rand(1,6);

% Valid constructor calls:
r = roi();
r = roi(x,y);
r = roi(x,y,s,p);
r = roi(x,y,[],[]);     % empty allowed
r = roi(x,y,s,p,img);
r = roi(x,y,s,p,img,ps,ipp,iop);

% Invalid constructor calls:
%roi(x)                          % X & Y must be provided together
%roi(x,y,s)                      % Slice & Phase must be provided together
%roi(x,y,s,p,img,ps)             %\_ PS, IPP, & IOP must be provided together
%roi(x,y,s,p,img,ps,ipp)         %/
assertExceptionThrown(@() roi(x),                  'ROI:roi:IncorectNargin')
assertExceptionThrown(@() roi(x,y,s),              'ROI:roi:IncorectNargin')
assertExceptionThrown(@() roi(x,y,s,p,img,ps),     'ROI:roi:IncorectNargin')
assertExceptionThrown(@() roi(x,y,s,p,img,ps,ipp), 'ROI:roi:IncorectNargin')


%% Test structured constructor method & version upgrading:
x = rand(10,1);
y = rand(10,1);
ps  = [0.5 0.5];        %\ 
ipp = [1 2 3];          % | These should be able to be 1xN or Nx1
iop = [1 0 0 1 1 1];    %/

% A version 0 structure:
s0 = struct(...
    'x',x,...
    'y',y,...
    'Slice',1,...
    'Interp','linear',...
    'Color','r',...
    'LineStyle','-',...
    'Tag','ROI-00');

% Version 1 structures:
s1a =  struct(...
    'x',x,...
    'y',y,...
    'Slice',1,...
    'DcmFile','/path/to/dicom/file.dcm',...
    'Interp','linear',...
    'Color','r',...
    'LineStyle','-',...
    'Tag','ROI-01');
s1b = s1a;
s1b.Phase = 1;

% Version 2 structure:
s2 = struct(...
    'x', x,...
    'y', y,...
    'Slice',1,...
    'Phase',1,...
    'Image','/path/to/image/file.dcm',...
    'Color','red',...
    'LineStyle',':',...
    'PixelSpacing', ps,...
    'ImagePositionPatient',ipp,...
    'ImageOrientationPatient',iop,...
    'Version', 2,...
    'Tag','ROI-02');

% upgrading Versions 0 and 1a should throw exceptions:
assertExceptionThrown(@() roi(s0),'ROI:struct2roi:DataRequired')
assertExceptionThrown(@() roi(s1a),'ROI:struct2roi:DataRequired')

% Old versions can be brought up to speed by minimally specifying 'Phase':
roi(s0,...
    'Phase',1);

% ...or perferably by specifying these three as well:
roi(s0,...
    'Phase',1,...
    'PixelSpacing',ps,...
    'ImagePositionPatient',ipp,...
    'ImageOrientationPatient',iop);

% And we can check we get errors if we incorrectly specify things when
% upgrading:
args = {s0,...
    'Phase',1,...
    'PixelSpacing',[0 0 0],...          % Incorrect PS - should be 1-by-2 or 2-by-1
    'ImagePositionPatient',ipp,...
    'ImageOrientationPatient',iop};
assertExceptionThrown(@() roi(args{:}),'ROI:checkpatientcs:IncorrectPixelSpacing')

args = {s0,...
    'Phase',1,...
    'PixelSpacing',ps,...
    'ImagePositionPatient',[1 2],...    % Incorrect IPP - should be 1-by-3 or 3-by-1
    'ImageOrientationPatient',iop};
assertExceptionThrown(@() roi(args{:}),'ROI:checkpatientcs:IncorrectImagePositionPatient')

args = {s0,...
    'Phase',1,...
    'PixelSpacing',ps,...
    'ImagePositionPatient',ipp,...
    'ImageOrientationPatient',[1 2 3]}; % Incorrect IOP - should be 1-by-6 or 6-by-1
assertExceptionThrown(@() roi(args{:}),'ROI:checkpatientcs:IncorrectImageOrientationPatient')

args = {s0,...
    'Phase',1,...
    'PixelSpacing',ps,...
    'ImagePositionPatient',ipp};        % Not enough inputs
assertExceptionThrown(@() roi(args{:}),'ROI:checkpatientcs:RequirePsIppIop')


% Upgrading version 1b should be elementary:
roi(s1b); 

% And N-dimensional structures should work too:
sNd(1) = s1b; sNd(2) = s1b;
r = roi(sNd);
assertEqual(size(sNd),size(r))



%% Test loading of all versions:
% Don't know what to do with version 0 or version 1a yet:
save('temp.mat','s0');
assertExceptionThrown(@() roi.loadfromfile('temp.mat'), 'ROI:struct2roi:DataRequired');
save('temp.mat','s1a');
assertExceptionThrown(@() roi.loadfromfile('temp.mat'), 'ROI:struct2roi:DataRequired');

% Version 1b and up is fine (they have "Phase"):
r = roi(s2); %#ok<NASGU>
save('temp.mat','s1b','s2','r')
[r,msg] = roi.loadfromfile('temp.mat');
assertTrue(isa(r,'roi'));   % Should have loaded correctly
assertEqual(3,numel(r));    % Should create 3 rois
assertTrue(~isempty(msg));  % But should return a message saying patient specs aren't defined


%% Test Tag conflicts:
r1 = roi(x,y);
r2 = roi(x,y);
save('temp.mat','r1')
R = load('temp.mat');
r3 = R.r1;
assertFalse(isequal(r1.Tag,r2.Tag)) % Two rois should have different tags
assertFalse(isequal(r1.Tag,r3.Tag)) % Saved & re-loaded rois get re-created tags


%% Test input arg size & shape

% Open / closed
r = roi(x,y);
r.isclosed;
r2 = r.close;
assertFalse(r.isclosed)
assertTrue(r2.isclosed)


%% Check Input size & shape tests
x = rand(1,10);
y = rand(1,10);

% Check that inputs get reshaped:
r = roi(x,y);
assertEqual(r.x,x(:))
assertEqual(r.y,y(:))

% Incorrectly sized arguments:
assertExceptionThrown(@() roi(x,x(1:end-1)),  'ROI:checkxy:InconsistentLengths')
assertExceptionThrown(@() roi(x(1:2),y(1:2)), 'ROI:checkxy:NotEnoughPoints')

% PS, IPP, & IOP
%              PixelSpacing  ->  1-by-2 float
%      ImagePositionPatient  ->  1-by-3 float
%   ImageOrientationPatient  ->  1-by-6 float
%
s = 1;
p = 1;

assertExceptionThrown(@() roi(x,y,s,p,img,ps,[1 2],iop),     'ROI:checkpatientcs:IncorrectImagePositionPatient')
assertExceptionThrown(@() roi(x,y,s,p,img,ps,[1 2 3 4],iop), 'ROI:checkpatientcs:IncorrectImagePositionPatient')

assertExceptionThrown(@() roi(x,y,s,p,img,ps,ipp,[1 2 3]),   'ROI:checkpatientcs:IncorrectImageOrientationPatient')
assertExceptionThrown(@() roi(x,y,s,p,img,ps,ipp,rand(1,7)), 'ROI:checkpatientcs:IncorrectImageOrientationPatient')

%% Check datatypes:

%assertExceptionThrown(@() roi(x,y,s,p,img,ps,false(3,1),iop),'ROI:checkpatientcs:IncorrectPixelSpacing')
%assertExceptionThrown(@() roi(x,y,s,p,img,ps,ipp,false(1,6)), 'ROI:checkixp:IncorrectImageOrientationPatient')


%% Test Dependent Getter methods:
x = [0 2 2 3 3 5 5 7 7 5 5 4 4 1 1 0 0];
y = [2 2 0 0 4 4 3 3 6 6 7 7 6 6 3 3 2];
r = roi(x,y);
L = 30;
A = 22;
% Correct values in pixels:
assertEqual(r.Length.pixels, L)
assertEqual(r.Area.pixels,   A)

% Conversion to mm not yet specified:
assertEqual(r.Length.mm, [])
assertEqual(r.Area.mm,   [])

% Provide coversion, and recalculate:
f = 0.5; 
r.PixelSpacing = [1 1]*f;
assertEqual(r.Length.mm, L*f)
assertEqual(r.Area.mm,   A*f^2)

%% Test some instance methods:

re = r.equispace(1);
figure, plot(r.x,r.y,'-*'), hold on, plot(re.x,re.y,'r*-')
keyboard



assertEqual( size(r.x), size(re.x) )
keyboard
dsv = hypot(diff(r.x),diff(r.y))



%% Test conversion to 3D
xy = [264.2465  258.0493  255.3760  254.5254  257.0772  264.8541  273.1171  277.9776  278.4637  276.5195  273.6031;
      187.8242  188.4318  191.8342  198.2745  203.6211  208.3601  207.9956  201.9199  194.7505  190.3760  187.3382];
xyz = [...
      -35.558769226073998  11.471424102784994  -6.200975799560595
      -35.558769226073998   8.445447540285002  -6.497655487060612
      -35.558769226073998   7.140125274659994  -8.158983612060609
      -35.558769226073998   6.724793243409991 -11.303661346435604
      -35.558769226073998   7.970789337159999 -13.914305877685607
      -35.558769226073998  11.768103790285011 -16.228270721435592
      -35.558769226073998  15.802771759034982 -16.050292205810607
      -35.558769226073998  18.176062774659982 -13.083641815185615
      -35.558769226073998  18.413416290284999  -9.582958221435604
      -35.558769226073998  17.464099884034994  -7.446971893310604
      -35.558769226073998  16.040076446534982  -5.963671112060609];
ps  = [  0.488281250000000   0.488281250000000 ];
ipp = [ -0.355587692260740  -1.170669059753400  0.850217781066894 ]*1e2;
iop = [  0 1 0 0 0 -1 ];


% Should fail if transformation parameters are not specified:
r = roi(xy(1,:),xy(2,:),1,1);
assertExceptionThrown(@() r.to3d, 'ROI:to3d:TransformationNotSpecified')

r.PixelSpacing = ps;
assertExceptionThrown(@() r.to3d, 'ROI:checkpatientcs:RequirePsIppIop')

r.ImagePositionPatient = ipp;
assertExceptionThrown(@() r.to3d, 'ROI:checkpatientcs:RequirePsIppIop')

r.ImageOrientationPatient = iop;

% But should pass and be correct if all transform parameters are specified
assertAlmostEqual(r.to3d, xyz)

% And preserve traces
X = r.to3dcell;
assertAlmostEqual( cat(1,X{:}), xyz )


%% Test joining / editing / inserting segments

% A basic job:
th = linspace(0,2*pi,100);
x1 = cos(th);
y1 = sin(th);
th = linspace(pi-pi/6,pi+pi/10,20);
x2 = cos(th)*2 + 1.5;
y2 = sin(th)*2;
figure, plot(x1,y1,x2,y2), axis equal
axis([-1 1 -1 1]*1.5)
r = roi(x1,y1);
r = r.close;
assertFalse(r.isclockwise)
assertTrue(r.toclockwise.isclockwise)
r = r.insert([x2(:) y2(:)]);
hold on
plot(r.x,r.y,'r-')

%% One that does not have either real or projected intersections:
s = 6;
np = 20*s;
th = linspace(pi/2-pi/6,pi/2+pi/6,np);
x1 = cos(th)*1.2*s;
y1 = sin(th)*s - s*.9;
r1 = roi(x1,y1);
r1 = r1.close;

x2 = linspace(-pi,pi,np);
y2 = cos(x2)+s*.2;
r2 = r1.join([x2(:) y2(:)]);

% One that has one real and one projected intersection:
x2 = linspace(-pi,pi/2,round(np*.75));
y2 = cos(x2)+s*.2;
r2 = r1.join([x2(:) y2(:)]);

% And then the other way around:
x2 = fliplr(x2);
y2 = fliplr(y2);
r2 = r1.join([x2(:) y2(:)]);

% And one that has two projected intersections:
x2 = linspace(-pi/2,pi/2,round(np*0.5));
y2 = cos(x2)+s*0.2;
r2 = r1.join([x2(:) y2(:)]);

% One that has two real intersections:
x2 = linspace(-pi,pi,np);
y2 = cos(x2)+s*0.15;
r2 = r1.join([x2(:) y2(:)]);



figure, plot(r1.x,r1.y,'b-+'), axis equal, hold on,
%plot(x2,y2,'r+')
plot(r2.x,r2.y,'r-+')




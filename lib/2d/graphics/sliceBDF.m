function sliceBDF(hObject,eventdata)

% Get the location of the click:
ax = gca;
xyz = surfacePoint( get(ax,'CurrentPoint'), get(hObject,'UserData') );

%plot3(gca,xyz(1),xyz(2),xyz(3),'r*')


% Get the relevant image:
handles = guidata(hObject);
s = find(round(handles.Images.z*1e3) == round(xyz(3)*1e3));
p = current('phase',handles.axes1);
img = handles.Images.image(s,p);


img = medfilt2(img);

% Get a Seed region:
figure, imshow(img,[]);
hg = imfreehand(gca);
mask = hg.createMask;
hg.delete;


% set input parameters
lambda = 0.5;
iterations = 600;

% perform segmentation
seg = sfm_chanvese(img,mask,iterations,lambda);
%seg = sparse(seg); % To save storage space

%I = imlincomb(0.5,single(img)/max(img(:)),0.5,single(seg));
%imshow(I,[])


[x,y] = mask2xy(seg);

hold on
%hf = fill(x,y,[0,1,0],'FaceAlpha',0.2,'Edgecolor','g');
contour(seg,[0 0],'r','linewidth',1)



% ------------------------------------------------------------------------
function [x,y,xfull,yfull] = mask2xy(mask)

mask = imdilate(mask,strel('diamond',1));
[r1,c1] = find(mask,1,'first');
B = bwtraceboundary(mask,[r1,c1],'NW');
xfull = B(:,2);
yfull = B(:,1);
% But now we need to simplify the result.
[x,y] = simplifyCurve(B);


% ------------------------------------------------------------------------
function [x,y] = simplifyCurve(C)

xf = C(:,2);
yf = C(:,1);

% extend the vectors by 50% each end
dn = floor(numel(xf)/2);
xe = [xf(end-dn+1:end) ; xf ; xf(1:dn)];
ye = [yf(end-dn+1:end) ; yf ; yf(1:dn)];
isorig = [false(dn,1); true(size(xf)); false(dn,1)];

xbase = xe(2:end-1);
ybase = ye(2:end-1);

% Compare point i to point i+1 and i-1
%   If the point is the same as both neighbours, it gets labelled as a
%   duplicate and removed.
isxdup = (xe(1:end-2)==xbase)  &  (xbase==xe(3:end));
isydup = (ye(1:end-2)==ybase)  &  (ybase==ye(3:end));

isdup = [false; (isxdup|isydup); false];

x = xe( ~isdup & isorig );
y = ye( ~isdup & isorig );


function [xyz] = surfacePoint(cpmatrix,zp)
% SURFACEPOINT Find (x,y,z) coordinates of a mouse-clicked point (cpmatrix)
%   with a constant z-surface with value zp.
%   The function calculates the (x,y,z) coordinates of the intersection of
%   the line described by cpmatrix and the point zp

dcp = diff(cpmatrix);               % (x,y,z) length of line
pct = (zp-cpmatrix(1,3)) / dcp(3);  % Percentage location of z between z1 & z2
xyz = cpmatrix(1,:)+dcp*pct;        % Cordinates of intersection



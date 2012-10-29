function newxy = axisloc2screenloc(ax,xy)
% AXISLOC2SCREENLOC Convert data locations on axis to Screen locations (pixels)
%
% ax - axis handle
% xy - N-by-2 set of points to convert
Transform = LocalAxis2Screen(ax);   
isLogX = strcmp(get(ax,'XScale'),'log');
isLogY = strcmp(get(ax,'YScale'),'log');
isRevX = strcmp(get(ax,'XDir'),'reverse');
isRevY = strcmp(get(ax,'YDir'),'reverse');
Xlim = get(ax,'Xlim');      % X limits
Ylim = get(ax,'Ylim');      % Y limits

% Compute normalized axis coordinates
X = xy(:,1);
Y = xy(:,2);

% Log-dependent calculations:
if isLogX,
    NormX = (log2(X) - log2(Xlim(1))) ./ (log2(Xlim(2))-log2(Xlim(1)));
else
    NormX = (X-Xlim(1)) ./ (Xlim(2)-Xlim(1));
end
if isLogY,
    NormY = (log2(Y)-log2(Ylim(1))) ./ (log2(Ylim(2))-log2(Ylim(1)));
else
    NormY = (Y-Ylim(1)) ./ (Ylim(2)-Ylim(1));
end

% Check for reversed axes:
if isRevX
    NormX = 1-NormX;
end
if isRevY
    NormY = 1-NormY;
end

XY = [NormX(:) NormY(:)];
r = size(XY,1);
A = repmat(Transform(1:2),r,1);
B = repmat(Transform(3:4),r,1);

newxy = A + B.*XY;
    
end  %axisloc2Screenloc()


% ------------------------------------------------------------------------
function T = LocalAxis2Screen(ax)
% Axis to screen coordinate transformation.
%
%   T = AXIS2SCREEN(AX) computes a coordinate transformation 
%       T = [xo yo rx ry] 
%   that relates the normalized axes coordinates [xa;ya] of a 
%   given point to its screen coordinate [xs;ys] (in the root 
%   units) by
%       xs = xo + rx * xa
%       ys = yo + ry * ya

% Get axes normalized position in figure
Fig = ancestor(ax,'figure');
AxisPos = hgconvertunits(Fig,get(ax,'position'),get(ax,'units'), ...
    'normalized',Fig);

% Get figure's normalized position in screen
FigPos = hgconvertunits(Fig,get(Fig,'position'),get(Fig,'units'), ...
    'normalized',get(Fig,'Parent'));


% Transformation norm. axis coord -> norm. fig. coord.
T = AxisPos;

% Transformation norm. axis coord -> norm. screen coord.
T(1:2) = FigPos(1:2) + FigPos(3:4) .* T(1:2);
T(3:4) = FigPos(3:4) .* T(3:4);

% Transform to screen units
ScreenSize = get(0,'ScreenSize');
T(1:2) = ScreenSize(1:2) + ScreenSize(3:4) .* T(1:2);
T(3:4) = ScreenSize(3:4) .* T(3:4);
end %LocalAxis2Screen()
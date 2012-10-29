function varargout = annotationManager(h,varargin)
%
% ANNOTATIONMANAGER(H,Value1,Value2,...) updates an existing annotation
% specified by the handle H.
%
% ANNOTATIONMANAGER(H) refreshes the position of the annotation H.
%
% H = ANNOTATIONMANAGER(HF,Position) creates a new empty annotation on the
%       figure specified by HF.
% H = ANNOTATIONMANAGER(HF,Position,Property1,Value1,...) overrides the
%       default properties.

switch get(h,'Type')
    case 'axes'
        
        if isequal(varargin{1},'NewCanvas')
            % annotationManager(axes,'NewCanvas',CanvasTag)
            hax = newCanvas(h,varargin{2});
            varargout{1} = hax;
        else
            han = newAnnotation(h,varargin{:});
            varargout{1} = han;
        end
        
        
    case 'text'  % custom annotation
        if isempty(varargin)        
            positionAnnotation(h)   % Update its position
        else
            updateAnnotation(h,varargin{:});
        end
        varargout = {};
        
    otherwise
        keyboard
        error('Unsupported object')
        
end 


% ------------------------------------------------------------------------
function hax = newCanvas(hax_ref,tag)

hp = get(hax_ref,'Parent');

hax = axes('Parent',hp,...
    'Units','normalized',...
    'Position',[0 0 1 1],...
    'Visible','off',...
    'HitTest','off',...
    'Tag',tag);


% ------------------------------------------------------------------------
function ht = newAnnotation(hax,tag,location,varargin)
assert(isequal(get(hax,'Type'),'axes'))
assert(any(strcmpi(location,{'NorthEast','NorthWest','SouthEast','SouthWest'})))

% Create the text object 
ht = text(0,0,'',...        % Origin, empty string for now
    'Parent',hax,...
    'Tag',tag,...
    'HitTest','off',...     % Don't capture mouse clicks
    'Interpreter','none',...
    'FontSize',11,...       
    'Color','w',...         % White font
    'LineStyle','-',...
    'EdgeColor','none',...
    'FontName','Helvetica',...
    'Visible','off',...
    'UserData',location);   % Store its location tag

% Position the annotation
positionAnnotation(ht); 

% Over-ride any extra properties:
if ~isempty(varargin)
    set(ht,varargin{:});
end


% ------------------------------------------------------------------------
function positionAnnotation(hAnno)
%POSITIONANNOTATION Position the annotation according to its location tag
%   which is stored in UserData.
% The possible location tags are:
%   'NorthEast','NorthWest','SouthEast','SouthWest'

%set(ht,'Units','pixels');

% Location test:
location = get(hAnno,'UserData');
lochas = @(str)~isempty(strfind(lower(location),str));

% Configure position:
set(hAnno,'Units','normalized');
buf = .01;                    % [pix] buffer


% Set horizontal position
if lochas('east')
    x = 1.0-buf;        % Locate in right half
    halign = 'right';   % Align left
elseif lochas('west')
    x = 0.0+buf;        % Locate in left half
    halign = 'left';    % Align right
end

% Set vertical position
if lochas('north')
    y = 1.0-buf;        % Locate in top half
    valign = 'top';
elseif lochas('south')
    y = 0.0+buf;        % Locate in bottom half
    valign = 'bottom';
end

set(hAnno,...
    'Position',[x,y,0],...
    'HorizontalAlignment',halign,...
    'VerticalAlignment',valign);


% ------------------------------------------------------------------------
function updateAnnotation(hAnno,varargin)
% 2D case modelled off Osirix, with some differences

switch get(hAnno,'Tag')
    case 'StackAnnotation'  % 3D annotation
        stackAnno(hAnno,varargin{:});
        
    case 'ImageAnnotation'  % 2D annotation
        imageAnno(hAnno,varargin{:});
       
    case 'PhaseText'        % 3D annotation
        phaseAnno(hAnno,varargin{:})
    
    otherwise
        error('Annotation not supported')
        
end


% ------------------------------------------------------------------------
function stackAnno(hAnno,zoom,curSlice,numSlices,curPhase,numPhases,Thickness,Location)
sanno = get(hAnno,'String');

% Initialisation case:
if isempty(sanno)
    sanno = {...
        'Zoom: N/A';
        sprintf('Slice: %2d/%2d',curSlice,numSlices);
        sprintf('Phase: %2d/%2d',curPhase,numPhases);
        sprintf('Thickness: %d mm  Location: N/A',Thickness)
        };
    
else
% Update case

    % Here we actually edit the existing string on a line-by-line basis, so
    % old info can remain by leaving those lines untouched, if desired.
    
    have = @(var)~isempty(var);
    
    % Update Zoom
    if have(zoom)
        sanno(1) = {sprintf('Zoom: %d%%',zoom)};
    end
    
    % Update Slice
    if have(curSlice) && have(numSlices)
        sanno(2) = {sprintf('Slice: %2d/%2d',curSlice,numSlices)};
    end
    
    % Update Phase
    if have(curPhase) && have(numPhases)
        sanno(3) = {sprintf('Phase: %2d/%2d',curPhase,numPhases)};
    end
    
    % Update Thickness & Location
    if have(Thickness) && have(Location)
        sanno(4) = {sprintf('Thickness: %5.2f mm  Location: %5.2f mm',Thickness,Location)};
    elseif have(Thickness)
        sanno(4) = {sprintf('Thickness: %5.2f mm',Thickness)};
    elseif have(Location)
        sanno(4) = {sprintf('Location: %5.2f mm',Location)};
    end
end

% Update string:
set(hAnno,'String',sanno);



% ------------------------------------------------------------------------
function imageAnno(hAnno,imsize,pt,value,PixelSpacing,ImagePositionPatient,ImageOrientationPatient)
set(hAnno,'Interpreter','tex');
ianno = get(hAnno,'String');

% Cursor value 
pt = pt-1;          % Dicom uses zero-based
pxyz = roi.pix2mm3D(pt,PixelSpacing,ImagePositionPatient,ImageOrientationPatient);
pt = round(pt);     % Whole mm values
imDimStr = @()sprintf('Image Size: %d \\times %d',imsize);
crsLocPix = @()sprintf('X: %d px  Y: %d px  Value: %0.0f',pt(1),pt(2),value);
crsLocMM  = @()sprintf('X: %04.2f mm  Y: %04.2f mm  Z: %04.2f mm',pxyz(1),pxyz(2),pxyz(3));

% Initialisation case:
if isempty(ianno)
    ianno = {...
        imDimStr();...       % Image Size: 420 
        '';...             % X: 154 px  Y: 205 px  Value: 145
        '';...             % X: 89.73 mm  Y: 124.87 mm  Z: 345.13 mm
        };
else
% Update case:

    % 
    have = @(var)~isempty(var);
    
    if have(imsize),  ianno(1) = {imDimStr()}; end
    
    if have(pt),
        
        % If cursor location is within the image bounds, show text
        %   If not, don't show text
        if ( (pt(1) >= 0) && (pt(1) <= imsize(1)) ) &&...
                ( (pt(2) >= 0) && (pt(2) <= imsize(2)) )
            ianno(2:3) = {crsLocPix(); crsLocMM()}; 
        else
            ianno(2:3) = {'',''};
        end
    end
end


% Update string
set(hAnno,'String',ianno)


% ------------------------------------------------------------------------
function phaseAnno(hAnno,j,n)
str = sprintf('Phase:  %2d / %d',j,n);
set(hAnno,'String',str,'interpreter','tex')

    
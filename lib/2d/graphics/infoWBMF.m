function infoWBMF(hf,~)
%INFOWBMF Window Button Motion Function which updates tha annotation info
handles = guidata(hf);

% Only makes sense to update if we have DICOM files loaded:
if isempty(handles.Images) || isempty(handles.Images.X)
    return
end

sj = current('slice',handles.axes1);
pj = current('phase',handles.axes1);

% Get cursor position:
cpsn = get(handles.axes1,'CurrentPoint');
cpsn = cpsn(1,1:2); % [x y]

% Get image dimensions:
%I = currentImage( handles, sj, pj ); <-- This is VERY SLOW
I = get(findobj(handles.axes1,'-depth',1,'Type','image'),'CData');


imsize = size(I);

% Get image value under cursor, limiting cursor to be in-bounds:
% Note that this value could be a single number (for grayscale image), or
% an rgb vector
loc = round(cpsn);
pmin = [ 1 1 ];
pmax = [ size(I,2) size(I,1) ];

loc = min([pmax; loc],[],1);
loc = max([pmin; loc],[],1);
v = I(loc(2), loc(1), :);
v = v(:);


% Scaling data:
ps  = handles.Images.info(sj,pj).PixelSpacing;
iop = handles.Images.info(sj,pj).ImageOrientationPatient;
ipp = handles.Images.info(sj,pj).ImagePositionPatient;


% Update the annotation:
annotationManager(handles.ImageAnnotation,...
    imsize,cpsn',v,ps,ipp,iop)


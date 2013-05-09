function infoWBMF(hf,~)
%INFOWBMF Window Button Motion Function which updates tha annotation info
handles = guidata(hf);

% Only makes sense to update if we have DICOM files loaded:
IMG = handles.( activeImageField(handles) );
if isempty(IMG) || isempty(IMG.X)
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

% Get image value under cursor, limiting to be in-bounds:
loc = round(cpsn);
pmin = [ 1 1 ];
pmax = [ size(I,2) size(I,1) ];

loc = min([pmax; loc],[],1);
loc = max([pmin; loc],[],1);
v = squeeze(I(loc(2), loc(1), :));


% DICOM-only data:
switch activeImageField(handles)
    case 'DICOM'
        ps  = handles.DICOM.info(sj,pj).PixelSpacing;
        iop = handles.DICOM.info(sj,pj).ImageOrientationPatient;
        ipp = handles.DICOM.info(sj,pj).ImagePositionPatient;
    case 'IMAGE'
        ps = NaN(2,1);   % Uness we manually put a calibration 
        ipp = NaN(1,3);
        iop = NaN(1,6);
end

% Update the annotation:
annotationManager(handles.ImageAnnotation,...
    imsize,cpsn',v,ps,ipp,iop)


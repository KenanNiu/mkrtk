function infoWBMF(hf,~)
%INFOWBMF Window Button Motion Function which updates tha annotation info
handles = guidata(hf);

% Only makes sense to update if we have DICOM files loaded:
if isempty(handles.DICOM) || isempty(handles.DICOM.X)
    return
end

sj = current('slice',handles.axes1);
pj = current('phase',handles.axes1);

% Get cursor position:
cpsn = get(handles.axes1,'CurrentPoint');
cpsn = cpsn(1,1:2);

% Get image dimensions:
imsize = size(handles.DICOM.X(:,:,sj,pj));

% Get image value under cursor:
v = [];
loc = round(cpsn);
if all(loc > [0 0]) && all(loc < imsize)
    v = handles.DICOM.X(loc(2), loc(1), sj, pj);
end

% Update the annotation:
annotationManager(handles.ImageAnnotation,...
    imsize,cpsn',v,...
    handles.DICOM.info(sj,pj).PixelSpacing,...
    handles.DICOM.info(sj,pj).ImagePositionPatient,...
    handles.DICOM.info(sj,pj).ImageOrientationPatient)


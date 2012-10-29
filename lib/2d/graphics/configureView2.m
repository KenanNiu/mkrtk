function configureView2(handles)
% CONFIGUREVIEW2 Configure the view.
%
% This function is called when a new set of DICOM files is loaded, or when
% a new session is loaded.
%
% See also UPDATESLICE

ha = handles.axes1;
s = current('slice',ha);
p = current('phase',ha);

if ( isempty(s) || isempty(p) ) ||...   % Image not yet displayed
        isempty(handles.DICOM)  ||...   % or image not loaded
        isempty(handles.DICOM.X) ||...  % or image empty
        s > size(handles.DICOM.X,3) || ...    % Or slice out of bounds
        p > size(handles.DICOM.X,4)           % Or phase out of bounds
    s = 1;
    p = 1;
end

% Display/hide graphics as required:
axTag = get(handles.axes1,'Tag');       % STORE
haveImage = ~( isempty(handles.DICOM) || isempty(handles.DICOM.X) );
if ~haveImage
    vis = 'off';
    I = [];
    clim = [0 1];
else
    vis = 'on';
    I = handles.DICOM.X(:,:,s,p);
    clim = handles.DICOM.CLim;
end

imshow(I,'Parent',ha);   % Show image
set(ha,'Tag',axTag)                             % Imshow cleared it
set(ha,'CLim',clim);
updateSlice(handles,s,p)
set(ha,'Visible',vis)
set([handles.StackAnnotation, handles.ImageAnnotation],'Visible',vis)


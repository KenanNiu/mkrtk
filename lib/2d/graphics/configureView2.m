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

IMG = handles.Images;       % Shorthand
[ns,np] = IMG.stacksize();  % Size of the stack

haveImage = ~( isempty(handles.Images) || isempty(handles.Images.X) );

% Check bounds on displayed slice & phase:
if ( isempty(s) || isempty(p) ) ||...   % Image not yet displayed
        ~haveImage || ...               % Or no image is loaded
        s > ns || ...                   % Or slice out of bounds
        p > np                          % Or phase out of bounds
    s = 1;
    p = 1;
end


% Default calibration tool state (used for non-DICOM images):
CALTOOL_STATE = 'off';

% Display/hide graphics as required:
axTag = get(handles.axes1,'Tag');       % STORE
if ~haveImage
    vis = 'off';
    I = [];
    clim = [0 1];
else
    vis  = 'on';
    clim = IMG.CLim;
    I    = IMG.image(s,p);
    
    % Account for the different configurations when reading image:
    switch IMG.Type
        case 'DICOM'
            
        case 'IMAGE'
            CALTOOL_STATE = 'on'; % Enable / show calibration tool
    end
    
end

imshow(I,'Parent',ha);   % Show image
set(ha,'Tag',axTag)      % Imshow cleared the tag
set(ha,'CLim',clim);
updateSlice(handles,s,p)
set(ha,'Visible',vis)
set([handles.StackAnnotation, handles.ImageAnnotation],'Visible',vis)
%set(handles.CalibrateImageTool,'Visible',CALTOOL_STATE,'Enable',CALTOOL_STATE)


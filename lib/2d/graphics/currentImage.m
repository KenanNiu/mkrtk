function [I,clim] = currentImage(handles,s,p)
%CURRENTIMAGE Get current image from the active image field
%
% HANDLES can have either of the following two image structures active:
%   handles.DICOM
%   handles.IMAGE
%
% This function uses the current slice and phase to correctly return the
% current image from the currently active image field.
%
% See also ACTIVEIMAGEFIELD

% Joshua Martin, 9-May-2013

if nargin < 3
    p = current('phase',handles.axes1);    
    if nargin < 2
        s = current('slice',handles.axes1);        
    end
    
end

fieldname = activeImageField(handles);
IMG = handles.(fieldname);

switch fieldname
    case 'DICOM'
        I = IMG.X(:,:,s,p);
    case 'IMAGE'
        I = IMG.X(:,:,:,s,p);
end

clim = IMG.CLim;

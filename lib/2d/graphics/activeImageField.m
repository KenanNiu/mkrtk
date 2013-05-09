function fieldname = activeImageField(handles)
%ACTIVEIMAGEFIELD return name of active image field
% 
% ACIVEIMAGEFIELD returns the name of the currently active image field in
% handles.  This will be either of the following:
%   'DICOM'
%   'IMAGE'
%
% In earlier versions of handles, the 'IMAGE' field does not exist, so this
% function is robust against that.
%
% NOTE: This function depends on the fact that the two structures
%       handles.DICOM and handles.IMAGE are used mutually exclusively.
%
% See also ACTIVEIMAGE

% Joshua Martin, 9-May-2013

if isfield(handles,'IMAGE') && ...
        ( ~isempty(handles.IMAGE.X) || ~isempty(handles.IMAGE.files) )
    fieldname = 'IMAGE';
elseif isfield(handles,'DICOM')
    fieldname = 'DICOM';
else
    fieldname = '';     % This should never be:
end

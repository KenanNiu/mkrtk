function [X,z,ps,dinfo] = readDicoms(varargin)
% READDICOMS Read DICOM files from cell array of full filepaths
%
% [X,z,ps,dinfo] = readDicoms( cellPthFile ) loads the dicom files listed
% in the N-by-1 cell array of full path-file strings
%
% This function is a wrapper for one of two methods for reading:
%   DCM4CHE     - Java 
%   DICOMREAD   - Matlab's native dicom reading utility

pathstr = [];
if iscell(varargin{1})
    % readDicoms( FullFileCellArray )
    files = varargin{1};
elseif ischar(varargin{1}) && iscell(varargin{2})
    % readDicoms( PathString, FileCellArray )
    [pathstr, files] = varargin{:};
else
    error('Incorrect inputs')
end

% Sort in natural order:
files = sort_nat(files);

% DCM4CHE is about 20-30% faster than DICOMREAD.  Hardly worth all the
% effort.   And it's getting the dicom header that is really slow.
% The full header read using the DCM4CHE class which uses java is actually
% about 4 times slower than using DICOMINFO, because DICOMINFO uses
% a mex file to do the basic parsing.  So the claim to fame here is that we
% don't bother reading in the full header - we only read in those fields
% that we are interested in.  The result is that the DCM4CHE method with
% minimal header retrieval is about 4 times faster than using MATLAB's
% DICOMREAD + DICOMINFO functions, but we've sacrificed reading much of the
% header to get that speed.
if dcm4che.libsloaded
    [X,z,ps,dinfo] = dcm4che_dicom_read(pathstr,files);
else
    %warning('readDicoms:dcm4cheNotFound',...
    %    'Could not find java files for dcm4che')
    [X,z,ps,dinfo] = native_dicom_read(pathstr,files);
end

end %readDicoms()


% ------------------------------------------------------------------------
function [X,z,ps,dinfo] = dcm4che_dicom_read(pathstr,files)

n = numel(files);
hw = waitbar(0,'Initialising...');
wstr = @(j)sprintf('Reading DICOM file %d of %d',j,n);

% We've already checked this, but just in case:
if ~dcm4che.libsloaded()
    error('Could not load dcm4che java libraries')
end

% Create a dcm4che object:
dcm = dcm4che([pathstr files{1}]);

% We can access the header fields individually, so instead of loading the
% whole header, let's keep it fast and just load the fields we need:
infofields = {...
    'SeriesDescription',...
    'SliceThickness',...
    'PixelSpacing',...
    'ImagePositionPatient',...
    'ImageOrientationPatient',...
    'NumberOfTemporalPositions',...
    'AcquisitionNumber',...
    'InstanceNumber'};
sargs = infofields;
[sargs{2,:}] = deal({});
dinfo = struct(sargs{:});


% Read first set of data:
X = dcm.getImage;
z = dcm.get('SliceLocation');
dinfo(1) = dcm.get(infofields);

% Expand for the remaining images:
X(:,:,2:n) = NaN;
z(2:n,1)   = NaN;
dinfo(2:n,1) = dinfo;

% Now read in the remaining images:
for j = 2:n
    waitbar(j/n,hw,wstr(j));
    % Load next file:
    dcm = dcm4che([pathstr files{j}]);
    dinfo(j) = dcm.get(infofields);
    X(:,:,j) = dcm.getImage();
    z(j) = dcm.get('SliceLocation');
end

ps = dcm.get('PixelSpacing');
close(hw)

end %dcm4che_dicom_read()


% ------------------------------------------------------------------------
function [X,z,ps,dinfo] = native_dicom_read(pathstr,files)
%NATIVE_DICOM_READ Read dicom images using Matlab's native DICOMREAD and DICOMINFO 
%

% Set up waitbar
n = numel(files);
hw = waitbar(0,'Gathering information...');
wstr = @(j)sprintf('Reading DICOM file %d of %d',j,n);

% Read in the first image & its specs:
dinfo = dicominfo([pathstr files{1}]);  % header info
ps = dinfo.PixelSpacing;                % 1-by-2, pixel spacing in mm
z = dinfo.SliceLocation;                % Record z-value
waitbar(0,hw,wstr(1));                  % Update waitbar
X = dicomread([pathstr files{1}]);      % Get the slice image


% Now expand variables to include all slices:
dinfo(2:n,1) = dinfo;       % Replicate as placeholders
z(2:n,1) = NaN;             % Fill with NaNs
X(:,:,2:n) = NaN;           % Fill with NaNs

% Read in the remaining slices:
for j = 2:n
    
    % Check that the next image is consistent with other images
    dinfo(j) = dicominfo([pathstr files{j}]);
    
    % Check stop criteria:
    %if ~isconsistent(dinfo(1),dinfo(j))
    %    warndlg(sprintf(...
    %        'Stack break: only %d of %d images read before finding inconsistent DICOM headers.',...
    %        j-1,n),'Stack Break','modal')
    %    error('Inconsistent DICOM headers')
    %end
    
    waitbar(j/n,hw,wstr(j));
    
    % Read the slice:
    X(:,:,j) = dicomread([pathstr files{j}]);
    z(j) = dinfo(j).SliceLocation;
    
end

close(hw)
end %native_dicom_read()




% ------------------------------------------------------------------------
function tf = isconsistent(dinfo1,dinfoj)
% Check to see if the two structures are part of the same scan:
fields = {'SeriesDescription','Width','Height','PixelSpacing'};

tf = false;
n = numel(fields);
for j = 1:n
    if ~isequal( dinfo1.(fields{j}), dinfoj.(fields{j}) )
        break
    end
end

if j == n
    tf = true;
end
end %isconsistent()

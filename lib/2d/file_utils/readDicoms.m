function [X,z,ps,dinfo] = readDicoms(varargin)
% READDICOMS Read DICOM files from cell array of full filepaths
%
% [X,z,ps,dinfo] = readDicoms( cellPthFile ) loads the dicom files listed
% in the N-by-1 cell array of full path-file strings
%
% This function is a wrapper for Matlab's native DICOMREAD function.
%
% It used to also interface with the custom dcm4ch3 class which used the
% dcm4che java libraries to read dicom files, but that was removed on
% 06-Mar-2012.  Explanation follows.
%
% Calling the DCM4CHE java libraries through the custom dcm4che matlab
% class object was about about 20-30% faster than DICOMREAD.  However,
% getting the dicom header was really slow, and we need that for spacial
% and temporal information. Reading the full header using the DCM4CHE class
% and java libs is actually about 4 times slower than using DICOMINFO,
% because DICOMINFO uses a mex file to do the basic parsing, and we have to
% do a lot of manual parsing of the data returned from java.  So the claim
% to fame here was that we could speed up the process and make it about 4
% times faster than the MATLAB DICOMREAD + DICOMINFO version by only
% reading in the dozen or so fields from that header that we were
% interested in.
%
% Considering all that, it isn't worth using the java libs if we sacrifice
% so much of the header info,  because in the future we may need more of
% the header information which would not be saved in older sessions.  So we
% now exclusively use the MATLAB functions again to load the entire header,
% then it's there if we save that Segmentation session and need to access
% that information later for upgrades or additional features.  Otherwise we
% would have to source the dicom files again and re-read their headers.


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


[X,z,ps,dinfo] = native_dicom_read(pathstr,files);


end %readDicoms()


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

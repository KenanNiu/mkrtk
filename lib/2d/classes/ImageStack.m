classdef ImageStack < handle
% Image stack class

%
% Property desctiptions:
%
%     pth: [string]      Path where images were loaded from
%   files: [cell]        1-by-S or S-by-P cell array of file names
%    Type: [string]      [ {DICOM} | IMAGE ] String identifying type of
%                        images loaded  
%    info: [struct]      Structure of image info. For Type='DICOM' it
%                        contains the info returned from DICOMINFO(), or
%                        for Type='IMAGE' it contains the following three
%                        fields for compatibility:
%                               PixelSpacing
%                               ImagePositionPatient
%                               ImageOrientationPatient
%       X: [NxMx(3)xSxP] Image data. For greyscale images (usually DICOM),
%                        X has size N-by-M-by-S-by-P; for truecolor images,
%                        X has size N-by-M-by-3-by-S-by-P
%       z: [SxP double]  Image slice locations (mm)
%       s: [2x1 double]  PixelSpacing
%    CLim: [1x2 double]  Image intensity limits
%

%TODO - Include load4Ddcm and loadImfiles in here?

properties
    pth = '';
    files = [];
    Type = '';      % 'DICOM' or 'IMAGE'
    info = [];
    X = [];
    z = [];
    s = [];
    CLim = [];
end

properties (Hidden=true, SetAccess=private)
    Version = 1;
end


methods
    
    function S = ImageStack(varargin)
        if nargin == 1 && isa(varargin{1},'struct')
            S = struct2class(varargin{1});
        end
        
    end %ImageStack()

    %----------------------------------------------
    function S = loadobj(S)
        % This method is called by default when any object of this class is
        % loaded.  If the automatic load succeeds, we have successfully
        % created the object.  If it fails, we will have a struct
        if isa(S,'struct')
            % Automatic load failed.
            % Need a try/ catch to inform the user that they might need to
            % upgrade with parameters?
            S = struct2class(S);
        end
    end %loadobj()
    
    %----------------------------------------------
    function I = image(S,slice,phase)
        % Get the selected image(s)
        switch S.Type 
            case 'DICOM'
                I = S.X(:,:,slice,phase);   % 4-D stack
            case 'IMAGE'
                I = S.X(:,:,:,slice,phase); % 5-D stack
        end
    end %image()
    
    %----------------------------------------------
    function varargout = stacksize(S)
        %STACKSIZE Return the number of slices & phases in the stack
        %   [ns,np] = stacksize(S)
        %     dims  = stacksize(S)
        %
       switch S.Type
           case 'DICOM'
               [~,~,ns,np] = size(S.X);
           case 'IMAGE'
               [~,~,~,ns,np] = size(S.X);
           otherwise
               ns = 0; np = 0;
       end
       if nargout < 2
           varargout = {[ns,np]};
       elseif nargout == 2
           varargout = {ns,np};
       else
           error('Incorrect nargout')
       end
    end %stacksize()
    
end %methods

methods (Static)
    
    %----------------------------------------------
    function S = LoadFiles(filelist)
        % Load Image files into an ImageStack object
        %
        %
        if ~iscell(filelist)
            filelist = {filelist};
        end
        S = eval(sprintf('%s()',mfilename)); % empty object
        
        if all( cellfun(@isimage,filelist) )               % Normal image formats selected
            S.Type = 'IMAGE';
            
            % Load normal image files
            [S.X,S.info] = loadImfiles(filelist);
            
            % Store other data:
            S.files = filelist;
            
            
        elseif all( cellfun(@isdicom,filelist) )           % Dicom image files selected
            S.Type = 'DICOM';
            
            % Load Dicom images:
            [S.X, S.z, S.s, S.info, S.files] = load4Ddcm(filelist);
        else
            error('One or more files selected had an un-recognised file type')
        end
        
        % Store shared properties:
        [pathname,~,~] = fileparts(filelist{1});    % Get default path
        S.pth   = pathname;                         % Store path
        S.CLim  = [ min(S.X(:)) max(S.X(:)) ];      % Store image limits
        
    end
    
    
end %methods (Static)

end %classdef

% ========================================================================
% Helper Functions


% ------------------------------------------------------------------------
function S = struct2class(s)
%STRUCT2CLASS Compatibility/converstion function
%
% This function can upgrade structs to the current version of the
% IMAGESTACK class.  
% Instigate empty class (should be a simple 1-by-1:
S = repmat( eval(sprintf('%s()',mfilename)) ,size(s));

% We'll assign only public or protected properties. Private properties
% should be set some other way, becasue there's a reason that they're
% private (like 'Version', which remains fixed)
%
% The following code block is a somewhat complex method of copying fields
% from the structure to the class object, but should be pretty
% comprehensive.  See also the ROI class for a very similar piece of code.

mc = metaclass(S);
if isprop(mc,'PropertyList')
    % 2011+ version, a convenient list is provided:
    names = {mc.PropertyList.Name};
    sa = {mc.PropertyList.SetAccess};
else
    % 2010 and earlier version, we must get the stuff manually
    names = getmcprop(mc.Properties,'Name');
    sa = getmcprop(mc.Properties,'SetAccess');
end
    %------------- Legacy helper function (for pre-2011)
    function list = getmcprop(props,field)
        list = cell(1,numel(props));
        for p = 1:numel(props)
            list{p} = props{p}.(field);
        end
    end %getmcprop()
    %-------------
pub  = strcmp(sa,'public');
prot = strcmp(sa,'protected');
names = names( pub | prot );

% Now process all elements of the structure and copy data across to the
% class object:
for k = 1:numel(s)
    tmp = upgradeStruct(s); % Check that we have a fully-fledged structure
    for j = 1:numel(names)
        S(k).(names{j}) = tmp.(names{j});
    end
end

end %struct2class()


% ------------------------------------------------------------------------
function s = upgradeStruct(s)
% Nothing to do yet
%
%
% VERSIONS:
%   0   - Structure, before the introduction of the class
%   1   - First edition of class. Class maintains all the same fields, but
%         adds a new one specifying the image stack type, either 'DICOM' or
%         'IMAGE'.  This comes contemporaneously with the ability to load
%         non-dicom images, so all previously save structures can be
%         directly upgraded to having this property set to 'DICOM'

% Version 0 -> Version 1
%   * add 'Type' field:
if ~isfield(s,'Version')
    s.Type = 'DICOM';
end

end %upgradeStruct()


% ========================================================================
% Image loading functions

% ------------------------------------------------------------------------
function [X,z,s,dinfo,Files] = load4Ddcm(Files)
% LOAD4DDCM Load DICOM images as a matrix with up to 4 dimensions.
%
% Inputs:
%   PATHNAME:  Cell array of full file path/name strings.
%
% Output:
%      X: N-by-M-by-S-by-P matrix where NxM is the image size, and S & P are the number
%           of slices and phases respectively
%      z: S-by-P matrix of SliceLocation header property
%      s: 2-by-1 PixelSpacing header property
%  dinfo: N-by-M-by-S-by-P matrix of structures containing some/all dicom
%           header properties
%  Files: S-by-P cell array of full file names in their appropriate
%           slice/phase positions.

if isa(Files,'cell')
    [X,z,s,dinfo] = readDicoms(Files);      % Read all files
    
else
    error('Incorrect inputs')
end

% Now reshape into a 4D stack where:
%   size(X) = [rows, cols, nslices, nphases]

% Work out how many different positions were scanned & how many slices,
% depending on what info we have available in the header (which is
% dependent on what scanner was used):

if isfield(dinfo,'NumberOfTemporalPositions') && ...
        isfinite(dinfo(1).NumberOfTemporalPositions);
    % PHILIPS defines this property, which is convenient:
    np = dinfo(1).NumberOfTemporalPositions;   % # positions
    if np == 0
        np = 1;
    end
    ns  = numel(dinfo)/np;                     % # slices in each position
    % Reshape:
    X = reshape(X,size(X,1),size(X,2),np,ns);
    X = permute(X,[1 2 4 3]);
    z = reshape(z,np,ns)';
    dinfo = reshape(dinfo,np,ns)';
    Files = reshape(Files,np,ns)';
    
elseif isfield(dinfo,'AcquisitionNumber')
    % SIEMENS seems to use AcquisitionNumber and InstanceNumber
    acqNum = [dinfo.AcquisitionNumber];
    % In this case, Acquisition Number stores the temporal position id of
    % each image
    if numel(unique(acqNum))>1
        % Sort:
        [~,ind] = sort([dinfo.InstanceNumber]);
        X = X(:,:,ind);
        z = z(ind);
        dinfo = dinfo(ind);
        % Get dim 3 & dim 4:
        np = numel(unique(acqNum));
        ns = sum(acqNum==acqNum(1));
        % Reshape - but in different order to above:
        X = reshape(X,size(X,1),size(X,2),ns,np);
        z = reshape(z,ns,np);
        dinfo = reshape(dinfo,ns,np);
        Files = reshape(Files,ns,np);
    end
    
end
end %load4Ddcm()


% ------------------------------------------------------------------------
function [X,iminfo] = loadImfiles(Files)
%LOADIMFILES Load normal image files as a matrix with up to 5 dimensions.
%
% Inputs:
%   PATHNAME:  Cell array of full file path/name strings, arranged in the
%           S-by-P order that you would like them displayed in. (see below)
%
% Output:
%        X: N-by-M-by-3-by-S-by-P matrix where NxM is the image size, 3 is
%           the number of image channels (usually 1 or 3) and S & P are the
%           number of slices and phases respectively 
%   IMINFO: S-by-P structure containing image information.  At this stage
%           it is filled to the appropriate size with NaNs.


if isa(Files,'cell')
    
    % Initialize variable(s)
    
    I = load_as_truecolor(Files{1});
    
    % Check that we're not going to load too much
    idata = whos('I');
    if idata.bytes/1024/1024 > 100 % ie, more than 100MB
        warning('You have selected to a large amount of data. Please be patient...')
    end
    
    % Stub out the full 5-D matrix:
    X = padarray( I, [0, 0, 0, size(Files)-1] );
    
    % Read all remaining files
    for j = 2:numel(Files)
        [s,p] = ind2sub(size(Files),j);
        X(:,:,:,s,p) = load_as_truecolor(Files{j});
    end
    
    % Fill out IMINFO structure with appropriately sized NaN arrays for now:
    iminfo = struct(...
        'PixelSpacing',             repmat( {NaN(2,1)}, size(Files)),...
        'ImagePositionPatient',     repmat( {NaN(1,3)}, size(Files)),...
        'ImageOrientationPatient',  repmat( {NaN(1,6)}, size(Files)));
    
else
    error('Incorrect inputs')
end


    %------------------------------------
    function I = load_as_truecolor(fname)
        % Load image as N-by-M-by-3, converting from indexed if necessary
        
        [cdata,cmap] = imread(fname);
        
        % Convert to true color, preserving class
        if ~isempty(cmap)
            type = class(cdata);
            
            % Switch to 1-based indexing for looking up colormap
            cdata = cdata+1;
            
            % Make sure cdata is in the range from 1 to size(cmmap,1)
            cdata = max(1,min(cdata,size(cmap,1)));
            
            % Extract r,g,b components
            f = double(intmax(type));
            r = zeros(size(cdata), type); r(:) = cmap(cdata,1)*f;
            g = zeros(size(cdata), type); g(:) = cmap(cdata,2)*f;
            b = zeros(size(cdata), type); b(:) = cmap(cdata,3)*f;
            
            I = cat(3,r,g,b);
            
        else
            I = cdata;
            
        end
    end %load_as_truecolor()

end %loadImfiles()


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
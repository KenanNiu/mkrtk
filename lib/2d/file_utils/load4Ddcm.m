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

% Now to re-shape....

% Work out how many different positions were scanned & how many slices:
ntp = dinfo(1).NumberOfTemporalPositions;   % # positions
ns  = numel(dinfo)/ntp;                     % # slices in each position

% This may need a more elegant/general approach:
% Reshape output to 4-D:
X = reshape(X,size(X,1),size(X,2),ntp,ns);
X = permute(X,[1 2 4 3]);
z = reshape(z,ntp,ns)';
dinfo = reshape(dinfo,ntp,ns)';
Files = reshape(Files,ntp,ns)';
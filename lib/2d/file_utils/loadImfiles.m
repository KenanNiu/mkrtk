function X = loadImfiles(Files)
%LOADIMFILES Load normal image files as a matrix with up to 5 dimensions.
%
% Inputs:
%   PATHNAME:  Cell array of full file path/name strings, arranged in the
%           S-by-P order that you would like them displayed in. (see below)
%
% Output:
%      X: N-by-M-by-3-by-S-by-P matrix where NxM is the image size, 3 is
%           the number of image channels (usually 1 or 3) and S & P are the
%           number of slices and phases respectively 

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
else
    error('Incorrect inputs')
end



% ------------------------------------------------------------------------
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

function I = iconmaker(imagefile)
% Helper function for stand-alone use in making icons for toolbars &
% pushbuttons

I = [];

% Get image data:
if exist('imagefile','var')==0
    [f,p] = uigetfile({'*.png;*.jpg;*.jpeg;*.gif;*.bmp;*.tif;*.tiff',...
        'Image Files (*.png,*.jpg,*.jpeg,*.gif,*.bmp,*.tif,*.tiff)'},...
        'Select an image/icon',fullfile(matlabroot,'toolbox','matlab','icons'));
    if isequal(f,0)
        disp('User cancelled')
        return
    end
    imagefile = [p,f];
    [X,map,alpha] = imread(imagefile);
    
elseif ischar(imagefile)
    [X,map,alpha] = imread(imagefile);
    
else 
    error('Unhandled')
end


imax = intmax(class(X));

if 1%all(alpha == imax) % no transparency set
    % Add iconEditor to path:
    addpath(fullfile(docroot,'techdoc','creating_guis','examples'))
    
    % Use iconEditor to do transparency:
    I = iconEditor('iconfile',imagefile);
    
else
    imax = double(imax);
    X = double(X);
    alpha = double(alpha);
    M = alpha==imax;
    M = repmat(M,[1,1,size(X,3)]);
    X(~M) = NaN;
    I = X./imax;
end


    
    



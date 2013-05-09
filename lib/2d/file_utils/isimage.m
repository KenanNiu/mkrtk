function tf = isimage(fname)
%ISIMAGE Test to see if file is probably an image 
%
% ISIMAGE calls IMFORMATS and inspects the file extension of the file
% specified to return true of false according to whether the file extension
% is registered as a legitimate image file format.

[~,~,ex] = fileparts(fname);
formats = imformats;
imexts = {formats.ext};
imexts = [imexts{:}];
imexts = cellsmash('.',imexts);
tf = any( strcmpi(imexts, ex) );

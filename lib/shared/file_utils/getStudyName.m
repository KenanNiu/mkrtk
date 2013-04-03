function studyname = getStudyName(pth)
%GETSTUDYNAME Get the study name from the filepath.  
%
% Input should be either:
%   a) a path containing a DICOM folder, or 
%   b) a path to a file used in registration
%
% Trailing file separator can either be present or not

c = path2cell(pth);



% Case sensitive: 'DICOM'
if any( strcmp(c,'DICOM') )
    idx = find( strcmp(c,'DICOM') );
    
    % Case insensitive: 'dicom'
elseif any( strcmpi(c,'dicom') )
    idx = find( strcmpi(c,'dicom') );
    
    % Case insensitive 'mgr'
elseif any( strcmpi(c,'mgr') )
    idx = find( strcmpi(c,'dicom') );
   
    % Case insensitive 'phantom'
elseif any( strcmpi(c,'phantom') )
    idx = find( strcmpi(c,'phantom') );
    
    % Case insensitive 'pilot'
elseif any( strcmpi(c,'pilot') )
    idx = find( strcmpi(c,'pilot') );
    
else
    
    % All known searches failed,
    %  return up to 3 levels of parent directories:
    if numel(c) >= 4
        idx = numel(c)-3:numel(c)-1;
    else
        idx = 1:numel(c)-1;
    end
    
end

% Insure it's a variable safe string
studyname = variablize(c(idx));

assert(~isempty(studyname),...
    'GETSTUDYNAME could not create a study name from the path')

function studyname = getStudyName(pth)
%GETSTUDYNAME Get the study name from the dicom filepath.  Input should be
% the path to the folder in which the DICOMDIR file resides.
%
% Trailing file separator can either be present or not

c = path2cell(pth);

% Take last segment of path:
studyname = c{end-1};
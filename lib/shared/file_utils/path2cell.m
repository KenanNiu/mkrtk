function pc = path2cell(pth)
%PATH2CELL Convert a filepath to a cell array of directories
%
% C = PATH2CELL(FILEPATH) converts FILEPATH to an 1-by-N cell array
% containing all the parent directories, with the final directory being
% C{end}.  If FILEPATH is a path to a file, the file does not appear in the
% cell array.
%
% PATH2CELL is platform independent, meaning that it will use the
% file separator that is present in FILEPATH, and not the one of the
% system.

% Find the file separator
if numel(strfind(pth,'\')) > numel(strfind(pth,'/'))
    delim = '\';
else
    delim = '/';
end

% Chop into segments:
if isequal(delim,'\')   % On Windows systems
    delim = '\\';       %  find '\' with '\\' using textscan
end
c = textscan(pth,'%s','Delimiter',delim);
pc = c{1}';

% Check for blank:
if isempty(pc{1})    % On unix systems we get an empty first entry
    pc(1) = [];      % Remove it
end




function make(mainfile)
% MAKE  Build executable for MKRTK
% USEAGE:  
%   MAKE
%
% ... not up to date with current build ...

if nargin == 0
    mainfile = 'MKRTK.m';
end

if ~ispc
    error('No ready for other platforms')
end


%===============================================================
warning('off','MATLAB:MKDIR:DirectoryExists')

%~~~~ Setup directory structure:
wkdir = pwd;    % Current directory
mkdir .Build    % Temp build folder
buildFolder = [wkdir filesep '.Build'];  % no trailing '\'

%~~~~ Binary folder:
cd .. % move up to program directory
standAloneDir = pwd;
cd(wkdir)

%~~~~ Supporting Files
%[supportfiles, iconfile, includefiles, fordelete] = dirstruct2string(dir(wkdir),mainfile);
mfiles = dirstruct2string(wkdir,mainfile,'*.m');
matfiles = dirstruct2string(wkdir,mainfile,'*.mat');
other = dirstruct2string(wkdir,mainfile,{'*.mex'});

%~~~~ Set up the build command.  Example:  
%   >> mcc -md 'C:\Documents and Settings\j.martin\My Documents\Stitch\Source\Build' stitch.m stitch.fig -M icon.res

compileString = ['mcc -md ''' buildFolder ''' ' mainfile ' '...
     mfiles matfiles other ];

%~~~~ Compile:
fprintf('Compiling file: %s . . . \n', mainfile)
eval(compileString)

%~~~~ Delete files converted for embedding:
%delete(fordelete{:})

% copy to "Standalone"
disp('===============')
try
    exe = [mainfile(1:end-2) '.exe'];  % Name of the compiled program
    copyString = ['!copy "', buildFolder, filesep, exe, '" "', standAloneDir, filesep, exe, '"'];
    eval(copyString)
    fprintf('File %s copied to:\n',exe)
    fprintf('<a href="matlab: eval([''! explorer %s''])">%s</a>\n',standAloneDir, standAloneDir)
    eval(['!rmdir /S /Q "' buildFolder '"'])
catch ME
    disp('There was an error when trying to copy to the Standalone directory:')
    rethrow(ME)
end

end %make()


%-------------------------------------------------------------------------
function files = dirstruct2string(root,mainfile,fspec)
% Create compiler string for supporting files from the directory structure
% DSTRUCT.

if ~iscell(fspec)
    fspec = {fspec};
end

if ~isequal(root(end),filesep)
    root(end+1) = filesep;
end

plist = textscan(genpath(root),'%s','Delimiter',pathsep);
plist = plist{1};

% Make relative:
nc = numel(root);
plist{1} = '';
for j = 2:numel(plist)
    plist{j}(1:nc) = [];
end

ishidden = @(p)~isempty(regexp(p,'\.*','ONCE'));
plist(cellfun(ishidden,plist)) = [];

files = {};
for j = 1:numel(plist)
    for e = 1:numel(fspec)
        if isempty(plist{j}) 
            pth = '';
        else
            pth = [plist{j} filesep];
        end 
        d = dir([pth fspec{e}]);
        % Collect and inject file separator:        
        F = cellfun(@(s)[pth s],{d.name},'UniformOutput',false);
        files = [files F];
    end
end

% Sort / unique
files = unique(files);
if isempty(files)
    files = ' ';
    return
end

% Drop the main file from the list:
files(~cellfun(@isempty,strfind(files,mainfile))) = [];

% Drop resource fork files:
files(~cellfun(@isempty,strfind(files,'\._'))) = [];

% Now convert to a string:
C(2,:) = files;
if isequal(fspec{e},'*.m')
    pre_s = ' ''';
else
    pre_s = ' -a ''';
end
C(1,:) = {pre_s};
C(3,:) = {''''};
C = C(:)';

files = cell2mat(C);

end %dirstruct2string()

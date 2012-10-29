function tf = download_fex_files(varargin)
%DOWNLOAD_FEX_FILES Download files from Matlab File Exchange
%
% DOWNLOAD_FEX_FILES(F1,F2,..) downloads the files specifed as string ids
%
% DOWNLOAD_FEX_FILES(C) downloads the files specified by a cell array of id
%                       strings 
%
% DOWNLOAD_FEX_FILES(S) downloads the files specified by the structure
%                       array S.  S is described below.
%
% DOWNLOAD_FEX_FILES(...,outputdir) places the files in the specified
%                       directory instead of the current.  The directory
%                       must already exist.
% The input structure S can take either of two forms.  The latest version
% of a file can be downloaded simply by specifying its id as a string:
%       S.id = '22846'
% or a specific version can be downloaded if the version and zip name are
% also specified:
%         S.id = '22846'
%        S.ver = '11'
%       S.file = 'matlab_xunit_3.1.zip'
%
% In the latter case, if the specified version cannot be found, an attempt
% will be made to get the current version and a warning issued.
%
% If a file cannot be downloaded, an error message will be printed.
%
% The output TF returns whether the file was downloaded and unziped
% correctely.
%
% Examples
%   download_fex_files({'22846','142317'}) % Download xUnit and findjobj
%
%
% Joshua Martin, 31 August 2012

[files,output_dir] = parse_inputs(varargin{:});

tf = download_helper(files,output_dir);

end %download_fex_files()


% ----------------------------------------
function [s,output_dir] = parse_inputs(varargin)
% Get the destination path:
if nargin > 1 && ischar(varargin{end}) && isdir(varargin{end})
    output_dir = varargin{end};
    varargin(end) = [];
else
    output_dir = [pwd filesep];
end
% Ensure format consistency:
if isequal(output_dir(end),filesep)
    output_dir(end) = [];
end

if isstruct(varargin{1})
    s = varargin{1};
else
    s = cell2struct(varargin,'id',1);
end
if ~isfield(s,'ver')
    [s.ver] = deal('');
end
if ~isfield(s,'file')
    [s.file] = deal('');
end
end %parse_inputs()


% ----------------------------------------
function tf = download_helper(files,target_dir)
%DOWNLOAD_FEX_FILES Download & extract files from FEX
%

nd  = numel(files);
tf  = false(nd,1);   % state of whether file was successfully retrieved
grv = true(nd,1);    % state of whether we got the requested version or not

% Now we use anonymous functions to programatically build the web
% addresses that we need.  Examples of these addresses are as
% follows.
% Simple download url:
%   http://www.mathworks.com/matlabcentral/fileexchange/37959?download=true
%
% If you specify a version, the version and filename must also be specified with the id:
%   http://www.mathworks.com/matlabcentral/fx_files/14317/19/findjobj.zip
%   
base_url = 'http://www.mathworks.com/matlabcentral/';
cur_url = @(s)[base_url 'fileexchange/' s.id '?download=true'];         % URL for latest file
ver_url = @(s)[base_url 'fx_files/' s.id '/' s.ver '/' s.file];     % URL for specific version

% Now try downloading and extracting all the files
%fw = max( cellfun(@numel,name) );
for j = 1:nd
    fj = files(j);
    ztemp = [target_dir filesep 'tmp.zip']; % Temporary zip file
    
    fprintf('Downloading file with id: %-*s ..... ',10,fj.id);
    GET_CURRENT = isempty(fj.ver);
    if ~GET_CURRENT
        [zfile,status] = urlwrite(ver_url(fj),ztemp);
        if ~status
            grv(j) = false;
            GET_CURRENT = true;
        end
    end
    if GET_CURRENT
        % Ok, we couldn't get the requested version, so let's try
        % getting the current version. We'll throw a little warning
        % as well, because we can't tell what the authors are going
        % to do with these programs in the future
        [zfile,status] = urlwrite(cur_url(fj),ztemp);
        if ~status
            fprintf('failed! %s.\n',fj.id)
            tf = false;
            % Skip this file & handle the problems in caller
            continue
        end
    end
    
    [~,fname,ext] = fileparts(zfile);
    iszip = isequal(ext,'.zip');
    if status 
        if iszip
            fprintf('unzipping ..... ')
            % Unzip the files into a directy by its id in its directory
            unzip(zfile,[ target_dir filesep fj.id filesep ]);
            fprintf('done!\n')
            tf(j) = true;
        else
            mvfile = [ target_dir filesep fname ext ];
            movefile(zfile,mvfile)
            fprintf(['\nWARNGING: The file %s was not a zip file and ',...
                'was not handled. Your function may not be installed.',mvfile])
        end
    end
    
end %for

% Clean up:
if exist(ztemp,'file') == 2
    delete(ztemp)
end

% Update paths:
addpath(genpath(target_dir));  % Add new folders which should have been created

% Shout out a warning if we couldn't connect to the net:
if all(tf==false)
    message_no_connection();
    
    % Or warn if we couldn't get specified versions
elseif any(grv==false)
    message_new_version(files(~grv));
    
end %if

    % ----------------------------------------
    function message_new_version(files)
        % Display a message if we had to download a new version
        fprintf('\n')
        fprintf('****************************************************************\n')
        fprintf('** WARNING: Check new versions for compatibility\n')
        fprintf('** \n')
        fprintf('** The FEX files listed below did not have the requested versions\n')
        fprintf('** available, so the latest versions were downloaded instead.\n')
        fprintf('** Note that these versions have not been checked for\n')
        fprintf('** compatibility.\n')
        for k = 1:numel(files)
            fprintf('**\t%*s   ==>   %s\n',8,files(k).id,ver_url(files(k)))
        end
        fprintf('****************************************************************\n')
        fprintf('\n')
    end %message_new_version()

    % ----------------------------------------
    function message_no_connection()
        % Display message that we didn't download any files
        fprintf('\n')
        fprintf('****************************************************************\n')
        fprintf('** UNABLE TO DOWNLOAD ANY OF THE REQUIRED FILES.\n')
        fprintf('** \n')
        fprintf('** Check your network connection settings and try again, or else\n')
        fprintf('** download the following files manually and place them in\n')
        fprintf('** the directory: \n**\t%s\n',target_dir);
        fprintf('** Files: \n')
        for j = 1:nd
            fprintf('**\t%s\n',cur_url(files(j)))
        end
        fprintf('****************************************************************\n')
        fprintf('\n')
    end %message_no_connection()


end %download_fex_files()

function check_dependencies(program)
% CHECK_DEPENDENCIES Check dependencies exist.
%   Check dependencies are on the file path.  If they are not, either
%       - download them (in the case of FEX functions)
%       - build them (if we have any functions that need mex-ing)
%       - OR:
%       - throw an error for mandatory functions
%       - throw a warning for desirable functions
%
%   Input argument PROGRAM flags the dependency check for :
%       - 2D segmentation   ==>  { '2d' | 'segmentation.m'   }
%       - 3D reconstruction ==>  { '3d' | 'reconstruction.m' }
%       

HARD = true;        % Program will break if we don't have it
SOFT = false;       % Loss of functionality, but will still work.



% --------- Mathworks File Exchange dependencies: ---------
segmentation_funs = {}; 
reconstruction_funs = { 'icp.m',         SOFT;   % For registration
                        'ICP_finite.m',  SOFT;   % For registration
                        'drawPlane3d.m', SOFT};  % For work in progress
shared_funs = { 'findjobj.m',   HARD;
                'getjframe.m', HARD;
                'sort_nat.m',  HARD };
% Run checks:
check_fex_deps(shared_funs)
switch lower(program)
    case {'segmentation.m', '2d'}
        check_fex_deps(segmentation_funs)
    case {'reconstruction.m', '3d'}
        check_fex_deps(reconstruction_funs)
end
        

% --------- Other dependencies: ---------
%   (none configured yet)

end

% ------------------------------------------------------------------------
function check_fex_deps(deps)
if isempty(deps)
    return
end
files = fex_id_lookup(deps(:,1));

% Expect FEX functions to exist here:
base_dir = 'lib/shared/externals/FEX Functions/';

n = size(files,1);

% Determine existence of files:
files_exist = @() 2 == cellfun(@exist,{files.mfile},repmat({'file'},[1,n]));
ex = files_exist();

if all(ex)
    return
end

% Download packages for any functions that don't exist:
ex(ex==0) = download_fex_files( files(~ex), base_dir );

% Assert that files of all HARD dependencies must exist
fail = [deps{:,2}] & ~ex;
if any( fail )
	% Failed.  Tell us which files we need to get manually:
    error_missing_files(fail)
end


    % ----------------------------------------
    function error_missing_files(tf)
        % Throw an error for the missing files flagged in TF
        base_url = 'http://www.mathworks.com/matlabcentral/';
        web_url = @(s)[base_url 'fileexchange/' s.id];
        str1 = sprintf('The following mandatory files could not be found nor downloaded:\n');
        str2 = sprintf('  %s\n',files(tf).mfile);
        str3 = sprintf('\nTry downloading them manually from the Matworks website:\n');
        str = {};
        for j = find(tf(:)')
            str{end+1} = sprintf('  %s',web_url(files(j))); %#ok<AGROW>
        end
        estr = [str1, str2, str3, str{:}];
        error(estr)
    end %error_missing_files()

end %check_fex_deps()


% ------------------------------------------------------------------------
function file_struct = fex_id_lookup(func_name)
% SOURCE_DATA Get source data for files on FEX.
%   
% To find the url of the zip files, download the zip then (on a Mac) do 
% "Get info" on the zip.
% Note that doing it this way the version number is fixed.  Mathworks
% doesn't keep all older versions, but they do keep some (they might keep
% the current one plus 10 previous versions or something? Or kept by date?)

file_data = {...
    'findjobj.m',   '14317','18', 'findjobj.zip'; 
    'getjframe.m',  '15830', '2', 'getjframe.zip';
    'icp.m',        '27804', '9', 'icp.zip';
    'ICP_finite.m', '24301', '1', 'icp_finite.zip';
    'sort_nat.m',   '10959', '4', 'sort_nat.zip';
    'drawPlane3d.m','24484','14', 'geom3d_2012.04.05.zip'};

[~,~,inds] = intersect(func_name(:,1),file_data(:,1));
file_struct = struct(...
    'id',file_data(inds,2),...
    'ver',file_data(inds,3),...
    'file',file_data(inds,4),...
    'mfile',file_data(inds,1) );
end %fex_source_data()

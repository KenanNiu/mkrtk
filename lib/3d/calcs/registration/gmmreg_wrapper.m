function transformed_model = gmmreg_wrapper(scene,model)
% Convenience wrapper for the different methods of calling GMMREG code
% which can be obtained here:
%   http://code.google.com/p/gmmreg/
%
% Two basic methods may exist:
%   (1) Compiled binary, eg, gmmreg.exe
%   (2) Matlab version with MEX file from mex_GaussTransform.c
%
% Option 1 (binary) is 4+ times faster and is more likely to return a
% finished result, so it is preferred if it is available. It does, however,
% require the existence of bone.ini for its configuration options. 
%
% Option 2 is avaliable for all platforms and will be compiled if a the
% source code exists but a platform specific mex file does not.


if gmmregBinary()
   transformed_model = gmmregBinary(scene,model);
   
elseif gmmregMex()
    transformed_model = gmmregMex(scene,model);
        
else
    error('There was a problem locating GMMREG.  Check your installation or try another solver')
end


end %gmmreg_wrapper


% ------------------------------------------------------------------------
function varargout = gmmregMex(scene,model)
% Usage:
% --------
%   tf = gmmregMex()                 <- Test for existence of files
%  xyz = gmmregBinary(scene,model)   <- Call matlab & mex code for registration

% Check if necessary files exist:
tf = ( 2==exist('gmmreg_L2','file') && ...         % Main function
    2==exist('gmmreg_L2_costfunc','file') && ...   % Necessary function in sub-directory
    (2==exist('mex_GaussTransform.c','file') || ...    % Require source file
    3==exist('mex_GaussTransform','file') ) ...     %  or compiled mex function
    ); % -Matlab implementation

% If that's all we're doing, return:
if nargin == 0
    varargout = {tf};
    return
end

% Check if the mex file exists, if not, build it:
if 3 ~= exist('mex_GaussTransform','file')
    build_GaussTransformMex;
end

% Setup and run:
config = initialize_config(model,scene,'rigid3d');
config.normalize = 1;
%[~,transformed_model,~,~] = gmmreg_L2(config);
[~,transformed_model] = my_gmmreg_L2(config);

varargout = {transformed_model};

end %gmmregMex()

% ------------------------------------------------------------------------
function varargout = gmmregBinary(scene,model)
% Usage:
% --------
%   tf = gmmregBinary()                 <- Test for existence of binary 
%  xyz = gmmregBinary(scene,model)      <- Call binary for registration
%

% Check if the binary exists:
if ispc
    ini  = 'bone.ini';          %\_ these two must be in the same directory AND on the path
    exec = 'gmmreg_demo.exe';   %/
    tf = (exist(exec,'file')==2) && (exist(ini,'file')==2);
else
    tf = false;
    % Don't have compiled Mac/Linux binaries yet
end

% If that's all we're doing, return:
if nargin == 0
    varargout = {tf};
    return
end

estr = sprintf('Executable file %s does not exist or is not on the MATLAB search path.',exec);
assert(tf,estr);

% We have the code, so prepare inputs:
disp('Writing temporary files...')
exec_dir = fileparts(which(exec));
opts = {'precision','%1.7e','delimiter','\t','newline','pc'};
dlmwrite([exec_dir filesep 'bone_model.txt'], model,opts{:});   % Write text file for model
dlmwrite([exec_dir filesep 'bone_target.txt'],scene,opts{:});   % Write text file for target

% The command must be run from the directory in which it resides
thisdir = pwd;
try
    disp('Running Gaussian Mixtures registration...')
    cd(exec_dir)
    cmd = sprintf('! %s %s %s ', exec, ini, 'rigid');
    eval(cmd)
    
    disp('Reading transformed model')
    transformed_model = dlmread('transformed_model.txt');
catch ME
    cd(thisdir)
    rethrow(ME)
end
cd(thisdir)

varargout = {transformed_model};

end %gmmregBinary()


% ------------------------------------------------------------------------
function build_GaussTransformMex
% Build MEX function from mex_GaussTransform.c

% Check the source c-file exists:
if 2 ~= exist('mex_GaussTransform.c','file')
    error('Could not find mex_GaussTransform.c for building');
end

% Switch to build directory:
wkdir = pwd;
cd(fileparts(which('mex_GaussTransform.c')))

% Build:
disp('Compiling MEX function: mex_GaussTransform.c ...')
mex mex_GaussTransform.c GaussTransform.c -output mex_GaussTransform
disp('Done!')

% Switch back to working directory:
cd(wkdir)

end %build_GaussTransformMex()


% ------------------------------------------------------------------------
function [param,transformed_model,history,config] = my_gmmreg_L2(config)
% Simplified and modified version of gmmreg_L2()

history.x = [ ];
history.fval = [ ];
if nargin<1
    error('Usage: gmmreg_L2(config)');
end
d = size(config.model,2); % number of points in model set
if (d~=2)&&(d~=3)
    error('The current program only deals with 2D or 3D point sets.');
end

x0 = config.init_param;
%
options = optimset('algorithm','active-set', 'display','iter', 'LargeScale','off','GradObj','on', 'TolFun',1e-010, 'TolX',1e-010, 'TolCon', 1e-10);
options = optimset(options, 'MaxFunEvals', config.max_iter);
options = optimset(options, 'GradObj', 'on');
param = fmincon(@gmmreg_L2_costfunc, x0, [ ],[ ],[ ],[ ], config.Lb, config.Ub, [ ], options, config);

transformed_model = transform_pointset(config.model, config.motion, param);
config.init_param = param;

end %my_gmmreg_L2()
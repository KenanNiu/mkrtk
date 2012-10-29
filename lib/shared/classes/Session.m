classdef Session < handle
%SESSION Abstract class for handling session files.
%
% SESSION is used for interacting with a GUI's handles data for saving and
% loading sessions.
%
% Currenty the supported invoking GUI's are:
%   Segmentation
%   Registration
%
% Methods:
%   Session.new     Create new session
%   Session.load    Load session from file
%   Session.save    Save session to file
%   Session.reset   Reset current session (restore to clean)
%   
% Usage:
%   HANDLES = Session.new(HANDLES) instantiates the required fields into
%   handles based on the invoking mfile (see above for supported GUIs)
%
%   HANDLES = Session.new(HANDLES,CALLER) explicitly specifies the method
%   to use.
%
%   HANDLES = Session.load(HANDLES,SFILE) loads the session *.mat file
%   SFILE and imports the appropriate fields into HANDLES based on the
%   invoking mfile.
%
%   HANDLES = Session.load(HANDLES,SFILE,CALLER) explicitly specifies the method
%   to use.
%
%   Session.reset('hard') re-launches the GUI of the invoking mfile after
%   warning the user.
%
%   HANDLES = Session.reset('soft',HANDLES) resets the appropriate fields
%   of HANDLES based on the invoking mfile.
%
%
% Joshua Martin, 12 Sep 2012

methods (Abstract)
    has_no_subclasses(obj) % We have to declare something...
end % methods (Abstract)


methods (Static)
    
    % ---------------------------------------------
    function handles = new(handles,opt)
        %NEW Create a new session by refreshing appropriate handles fields
        %
        % If OPT is not specified, defaults to invoking filename
        if ~exist('opt','var')
            ds = dbstack(1);
            [~,opt,~] = fileparts(ds(1).file);
        end
        dflt_cell = defs(opt,'all');
        hs = cell2struct(dflt_cell(:,2),dflt_cell(:,1));
        handles = load_fields(handles,hs,fieldnames(hs),0,1);
    end %new()
    
    % ---------------------------------------------
    function [handles,ok] = load(handles,sfilename,opt)
        %LOAD Helper function to load session files
        %
        % LOAD(HANDLES,SFILE) loads the specified session file SFILE (*.mat
        % file containing handles data) into HANDLES and returns the
        % updated HANDLES structure.  The fields loaded in from SFILE are
        % predfined in this function.  Which set of fields to load is
        % selected based on the name of the calling mfile.
        %
        % LOAD(HANDLES,SFILE,OPT) explicitly sets which option to use for
        % loading.  See the case statement in the body of the function for
        % options.
        
        ok = false;
        
        % Load data
        if ~exist(sfilename,'file')
            errordlg('File does not exist')
            return
        end
        container = load(sfilename);
        
        % Check validity - must have field HANDLES:
        if ~isfield(container,'handles')
            errordlg('This file is not a session file.')
            return
        end
        
        % Get saved handles data:
        h = container.handles;
        
        % If OPT is not specified, try using the invoking mfilename
        if ~exist('opt','var')
            ds = dbstack(1);
            [~,opt,~] = fileparts(ds(1).file);
        end
        
        handles = load_fields(handles,h,defs(opt,'commonfields'),0);
        [handles,count] = load_fields(handles,h,defs(opt,'uniquefields'),0);
        
        % Now warn if we had problems:
        if count == 0
            warndlg('It appears that this was not a valid session file for this program')
            return
        end
        [p,~,~] = fileparts(sfilename);
        handles.sessionPath = p;
        ok = true;
                
    end %load()
    
    
    % ---------------------------------------------
    function varargout = reset(method,varargin)
        %RESET
        %
        % RESET('hard')
        % handles = RESET('soft',handles)
        varargout = {[]};
        if ~exist('method','var')
            method = 'hard';
        end
        % Warn user:
        btn = warncanceldlg('This will remove all data from your current session. Are you sure you want to continue?','Erase session?');
        if ~isequal(lower(btn),'ok')
            return
        end
        ds = dbstack(1);
        [~,caller,~] = fileparts(ds(1).file);
        % Choose mode of reset:
        switch lower(method)
            case 'hard'
                % For a hard reset, we delete the callback gui, then
                % re-instantiate it using the calling file.
                delete(gcbf);
                feval(caller);
            case 'soft'
                % For a soft reset, we get handles, re-initialise the
                % fields, then push back to gui.
                handles = varargin{1};
                handles = Session.new(handles,caller);
                guidata(gcbf,handles);
                varargout = {handles};
            otherwise
                error('not yet implemented')
        end
    end %reset
    
    % ---------------------------------------------
    function [filepath,ok] = save(defaultpath,handles) %#ok<INUSD>
        %SAVE Shared helper function to save session files
        %
        % SAVE saves all data in HANDLES, using DEFAULTPATH as the default
        % path for the UIPUTFILE dialog.
        
        % Request mat file name
        [filename, pathname] = uiputfile( ...
            {'*.mat','MAT-files (*.mat)'},'Save session as',defaultpath);
        
        % Abort if the user cancelled:
        ok = false;
        if isequal(filename, 0)
            filepath = '';
            disp('User Cancelled')
            return
        end
        
        filepath = [pathname filename];
        
        % Not sure why, but sometimes the save will crack the sads when
        % trying to over-write an existing file.  If this is the case, give
        % it a couple of shots with a pause in between to see if it will
        % work eventually:
        n = 5;      % Try up to 5 times
        for j = 1:n
            try
                disp('Saving session...')
                save(filepath,'handles')
                disp('Session saved successfully.')
                ok = true;
                break
            catch ME %#ok<NASGU>
                pause(0.5)
                if j == n
                    errordlg('There was an error and your data was not saved')
                else
                    disp('Failed... trying again...')
                end
            end
        end %for
    end %save()

end % methods (Static)


end % classdef


% ========================================================================
function out = defs(opt,type)
%
% This is a horrible method of storing defaults.
%
%  OPT = [ 'segmentation' | 'registration' ]
% TYPE = [ 'fields' | 'uniquefields' | ]

COMMON = {...
    'sessionPath', '';...
    'userPath',    [pwd filesep];...
    };

if strcmpi(type,'commonfields')
    out = COMMON(:,1);
    return
end

switch lower(opt)
    case 'segmentation'
        D.pth   = [];    % Initialise
        D.files = [];
        D.info  = [];
        D.X     = [];
        D.z     = [];
        D.s     = [];
        D.CLim  = [];   % Display limits of X
        UNIQUE = {...
            'DICOM',  D;...
            'traces', roi([]);...   % Empty roi structure
            };
    case {'registration','reconstruction'}
        UNIQUE = {...
            'Models',       []; ...
            'motion',       [];...
            'HelicalAxis',  [];...
            'LineOfAction', [];...
            'MomentArm',    [];...
            };
end %switch

if strcmpi(type,'uniquefields')
    out  = UNIQUE(:,1);
    return
end

out = cat(1,COMMON,UNIQUE);
if strcmpi(type,'fields')
    out = out(:,1);
end

end

% ========================================================================


% ------------------------------------------------------------------------
function [handles,count] = load_fields(handles,source,updfields,warn,force)
% Load the specified fields into handles, checking for existence in both
if ~exist('warn','var')
    warn = 1;
end
if ~exist('force','var')
    force = 0;
end
count = 0;
for j = 1:numel(updfields)
    if isfield(source,updfields{j})         % is field in source
        if isfield(handles,updfields{j})... % is field in destination
                || force                    % or force it
            handles.(updfields{j}) = source.(updfields{j});
            count = count+1;
        end
    else
        if warn
            fprintf('Warning: HANDLES field not found in session file: %s\n', updfields{j});
        end
    end
end
% Versioning issues - this is a bit ugly and could be dropped in due course
handles = check_cloud_upgrade(handles);
handles = check_roi_upgrade(handles);
end %load_fields

% ------------------------------------------------------------------------
function handles = check_cloud_upgrade(handles)
% Clouds used to be structs, so here we upgrade old structs to classes:
if isfield(handles,'Models') && isequal(class(handles.Models),'struct')
    M = handles.Models;
    % "HiRes" & "LoRes" clouds might need to be updated.
    % It would be more ambiguous do this in vectorised form across all
    % models because there may be some that don't exist, so we would have
    % to test for that, but then convert empty structures to empty objects.
    % Easier to just loop:
    for mj = 1:numel(M)
        if isequal(class(M(mj).HiRes),'struct')
            M(mj).HiRes = Cloud(M(mj).HiRes);       % Convert HiRes
        end
        if isequal(class(M(mj).LoRes),'struct')
            M(mj).LoRes = Cloud(M(mj).LoRes);       % Convert LoRes
        end
    end
    % "Zeroed" clouds no longer used
    if isfield(M,'Zeroed')
        % And because they aren't used anymore, we also need to either
        % trash or convert their pose info, if it's not empty.  We'll be
        % nice and convert it...
        
        % First, convert to Cloud:
        % (Zeroed clouds were dropped before implementing the class, so
        % they will always be structs)
        for mj = 1:numel(M)
            if ~isempty(M(mj).HiRes)
                M(mj).Zeroed = Cloud(M(mj).Zeroed);     % Convert Zeroed to clouds
                if isfield(M(mj),'q')
                    zp = M(mj).Zeroed.transform(M(mj).q.inverse,M(mj).x);   % Positioned clouds
                    sr = M(mj).HiRes(ones(size(zp)));               % Static, replicated
                    [M(mj).q,M(mj).x] = gettransform(sr,zp,'quat'); % new transformation
                end
                if isfield(M(mj),'qraw')
                    % Same as above, just with the 'raw' fields:
                    zp = M(mj).Zeroed.transform(M(mj).qraw.inverse,M(mj).xraw);   % Positioned clouds
                    %sr = M(mj).HiRes(ones(size(zp)));               % already exists above
                    [M(mj).qraw,M(mj).xraw] = gettransform(sr,zp,'quat'); % new transformation
                end
            end
        end
        % Ok, done.  Messy, but not too painful. Now delete zeroed clouds:
        M = rmfield(M,'Zeroed');
    end
    handles.Models = M;
end
end

% ------------------------------------------------------------------------
function handles = check_roi_upgrade(handles)
% Here we have to handle some versioning problems with rois:
if isfield(handles,'traces') && isa(handles.traces,'struct')
    % Upgrade traces:
    if isempty(handles.traces)
        handles.traces = roi({});
    else
        traces = roi(handles.traces);
        handles.traces = traces.addpatientcs(handles.DICOM.info);
    end
end
end %check_roi_upgrade()


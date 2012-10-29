function [qout,xout,flag] = registerSequence(models,q,x,solver,progfun)
% We ignore the q & x fields which are in Models
%
%

% Initialize:
qout = q;
xout = x;
flag = 0; 

% Abstract the clouds away from the models structure.  We must use cells
% because some of these will not exist and so we just skip them
static = {models.HiRes};
dynamic = {models.LoRes};

% Number of phases / objects:
np = max( cellfun(@numel,{models(:).LoRes}) );  % Number of phases
no = numel(models); % Number of objects

% Work out if we have clouds required for calculations:
haveClouds = ~cellfun(@isempty,static).*~cellfun(@isempty,dynamic);
haveClouds = ones(np,1)*haveClouds; % matrix form

% Now create cleverly ordered vectors containing the indices of the
% phase/object combinations that we need to analyse:
[plist,olist] = makeList(q.isnan .* haveClouds);

% Total number of registrations to perform:
nr = numel(plist);  % Total number of registrations to perform
rj = 0;             % Current registration

% Set up waitbar if no progress function has been supplied
wstr = @(rj,pj,oj)sprintf('Operation %d of %d:  Bone %d/%d;  Position %d/%d',rj,nr,oj,no,pj,np);
if ~exist('progfun','var') || ~isa(progfun,'function_handle')
    hw = specialWaitbar(wstr(0,0,0));
    progfun = @(x,s)waitbar(x,hw,s);    
end

% Set up our termination flag:
%   ICPSTOP == 0    => Continue calculations
%   ICPSTOP == 1    => Terminate and save current calcs
%   ICPSTOP == -1   => Terminate and discard current calcs
global ICPSTOP
ICPSTOP = 0;

% Calculate the transformations of the static clouds into the position of
% the dynamic clouds for every phase
t = tic;
for j = 1:nr
    
    pj = plist(j);
    oj = olist(j);
    
    % Check if processing is required on this object/phase:
    if ~q(pj,oj).isnan
        % Object already solved in this phase - skip it:
        continue
    end
    
    % Check for termination:
    pause(eps)  % Force matlab to process other threads (ie, our java button callbacks in progfun)
    if ICPSTOP == 1
        % "Stop" button pressed
        % Return calcs up to here
        break
    elseif ICPSTOP == -1
        % Waitbar's "Cancel" button pressed OR Waitbar destroyed:
        % Return the inital (q,x) unchanged
        return
    end
        
    % Now calculate the mapping from the static clouds into the
    % position of the dynamic clouds.  Q and X will then be used by
    % applying them to static clouds to find the reconstructed position.
    %
    % The procedure for doing this, however, is not to directly use the
    % static cloud in the matching process.  We instead use an initial
    % guess for the position, to give the matching procedure a better
    % chance at getting the correlation correct. 
    
    % Get point sets:
    static0 = static{oj};
    target  = dynamic{oj}(pj);
    
    % Get & apply an initial pose estimate to the hi-res model:
    [qi,xi] = getPoseEstimate;
        
    % And get pose:
    [qj,xj] = registerClouds(solver,static0,target,qi,xi);
    q(pj,oj) = qj;
    x(pj,oj,:) = xj;
    
    % Increment registration counter (for waitbar & overall progress)
    rj = rj+1;
    
    % Update waitbar.
    if toc(t) > 0.1 %nr < 100 || ( nr < 1000 && mod(j,10)==0 ) || mod(j,100)==0
        t = tic;
        %waitbar((rj-1)/nr,hw,wstr(rj,pj,oj))
        progfun( (rj)/nr, wstr(rj,pj,oj) );
    end
    
end %for

if exist('hw','var') && ishandle(hw)
    close(hw)
end

qout = q;
xout = x;
flag = 1;



    % -----------------------------------------------
    function [q0,x0] = getPoseEstimate
        % Get an estimate for the 
        tmp = double(~q(:,oj).isnan);
        tmp(tmp==1) = find(tmp==1);   % convert 1's to index
        tmp(tmp==0) = Inf;            % convert 0's to Inf
        
        % Now find the nearest measurement:
        [dist,idx] = min(abs(tmp-pj));      % nearest measurement
        
        % If found, use the nearest as seed:
        if dist ~= Inf
            q0 = q(idx,oj);
            x0 = squeeze(x(idx,oj,:));
        else
            % If no nearby pose was found, retrieve from static models:
            warning('registerClouds:NoSeed',...
                ['Registration is starting with unsupervised initial position.',...
                ' Be sure to check solution.'])
            
            q0 = quaternion([1 0 0 0]);
            x0 = [0 0 0];
        end
    end


end


% ------------------------------------------------------------------------
function [plist,olist] = makeList(qisnan)
%MAKELIST setup phase list and object list for processing
%
% The phase list PLIST and object list OLIST are vectors of corresponding
% phase/object indices that need to be analysed.  The analysis loop should
% work progressively through these vectors.
%
% The simple solution here would be to create the lists based off indices
% from the q.isnan matrix:
[plist,olist] = find(qisnan);
%
% However, we would like to order each object list so that we start with
% the first object that has been analysed, then work through to the end,
% then return to this point and work backwards.  So essentially, work
% outwards from the first analysed location.  This way the solving should
% be more robust, and we can initialise on any object/phase combination we
% like.
%
plist(:) = NaN;
[~,rstart] = min(qisnan,[],1);   % Find first analysed row in each col.

% Re-order each column separately:
r1 = 0;
for j = 1:size(qisnan,2)
    plistj = find(qisnan(:,j)); % default sorting
    
    % Now re-order by truncating:
    plistj = [plistj(plistj>=rstart(j)); flipud(plistj(plistj<rstart(j)))];
    
    % Update indices:
    r0 = r1+1;
    r1 = numel(plistj)+r0-1;
    
    % Now populate the vectors with re-ordered indices:
    plist(r0:r1) = plistj(:);
    olist(r0:r1) = j;
        
end

end


% ------------------------------------------------------------------------
function hw = specialWaitbar(str)
%SPECIALWAITBAR  Create a custom waitbar which allows stopping & cancelling
%   using global variables.  
% 

hw = waitbar(0,str,'CreateCancelBtn',@cancelfun,'Visible','off');
psn = get(hw,'Position');
h = 17;
w = 50;
s = 10;
v = 6;
delete(findobj(hw,'Tag','TMWWaitbarCancelButton'))
uicontrol(hw,'Style','pushbutton',...
    'String','Stop',...
    'TooltipString','Stop processing here & return results',...
    'Units',get(hw,'Units'),...    % Same as waitbar
    'Position',[(psn(3))/2-w-s/2, v, w, h],...
    'Callback',@stopfun)

uicontrol(hw,'Style','pushbutton',...
    'String','Cancel',...
    'TooltipString','Cancel processing & discard results',...
    'Units',get(hw,'Units'),...    % Same as waitbar
    'Position',[(psn(3))/2+s/2, v, w, h],...
    'Callback',@cancelfun)
set(hw,'Visible','on')
end %specialWaitbar

% ------------------------------------------------------------------------
function stopfun(~,~)
global ICPSTOP;
ICPSTOP=1;
delete(gcbf);
%drawnow
end %stopfun

% ------------------------------------------------------------------------
function cancelfun(~,~)
global ICPSTOP;
ICPSTOP=-1;
delete(gcbf);
%drawnow
end %cancelfun

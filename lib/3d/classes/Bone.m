classdef Bone
% BONE class.  
%
% A BONE class object represents a bone model and can store its name
% (.Tag), its high-resolution representation derived from a static 3D scan
% (.HiRes), its low-resolution representation derived from a dynamic scan
% (.LoRes), and the pose states which describe the mapping of the
% high-resolution model into the positions of the dynamic low-resolution
% model in terms of its rotation quaterion (.q) and translation vector
% (.x).  
%
% The rotation and translation (q and x) are outputs of the registration
% algoritm, and the BONE object can store both the raw result, and a
% smoothed result.  These raw/smoothed states are internally managed so
% that retrieving them returns the smoothed state if it exists, or the raw
% state if it does not.  The smoothed and raw states are stored in hidden
% properties with private set access (see below) so they can be retrieved
% if necessary, but normal setting/getting of 
%
%
%
% - Internally manages the raw/smoothed states of q & x
% - Supersedes "Models" structure.
%
% Properties:
%   Tag         (Ordinary) Bone tag or name
%   HiRes       (Ordinary) High Resolution (static) cloud objects
%   LoRes       (Ordinary) Low Resolution (dynamic) cloud objects
%   smoothed    (Ordinary) True/false - whether q & x return smoothed data
%   q           (Dependent) Returns qsmooth if it is not empty, otherwise qraw
%   x           (Dependent) Returns xsmooth if it is not empty, otherwise xraw
%
% Hidden Properties:
%   qraw        Quaternion orientation states as calculated from registration
%   xraw        Position states ([x y z]) as calculated from registration
%   qsmooth     Smoothed quaternion states (if smoothing has been performed)
%   xsmooth     Smoothed position states (if smoothing has been performed)
%   Version     Class version number (for backward compatibility)
%   
%
% Constructor:
%   b = Bone()                  % Create an 1-by-1 blank bone object%
%
% State smoothing methods:
%   b = b.smoothpose(fun)       % Smooth pose states (qraw & xraw) with the
%                               % function FUN, where FUN is a function 
%                               % handle
%   b = b.clearsmoothing        % Remove the smoothing so that the   
%                               % .q and .x methods will return the states 
%                               % .qraw and .xraw
%
%
% Because "q" and "x" are dependent properties, they have associated getter
% and setter methods.  They can therefore not be vectorised and should be
% called on arrays in the following manner:
%
%   Q = [b.q]                       % Getter methods
%   X = [b.x]                       %
%
%   [b.q] = deal(q1,q2,q3,...)      % Setter methods
%   [b.x] = deal(x1,x2,x3,...)      %
%
% The getter methods retrieve the smoothed states, qsmooth & xsmooth, if
% they are not empty, otherwise they return the raw states, qraw & xraw.
% 
% The setter methods for q and x reset the raw states, qraw & xraw, to the
% values specified and reset the smoothed states to empty.  Setting of the
% smoothed/raw states directly is currently not permitted since the
% smoothing should be performed by using the state smoothing methods
% provided.
%
% See also CLOUD

% Joshua Martin, 22-Nov-2012

properties
    Tag = ''
    HiRes
    LoRes
    smoothed = false
end

properties (Dependent=true)
    q
    x
end

properties (Hidden=true, SetAccess=private)
    qraw
    xraw
    qsmooth
    xsmooth
    Version = 1
end


methods
    
    % ------------------------------------------
    function b = Bone(varargin)
        % Constructor:
        %   Bone()          1-by-1 empty Bone object
        %   Bone(struct)    Bone objects converted from structs
        if nargin == 0
            % Create blank Bone object, ok
            return
            
        elseif nargin == 1 && isa(varargin{1},'struct')
            % Convert from struct to bone object
            b = struct2obj(varargin{1});
            
        else
            error('Incorrect input arguments')
            
        end
    end %Bone()
    
    
    % ------------------------------------------
    function b = clearsmoothing(b)
        % Remove smoothing and revert to raw pose states
        [b.smoothed] = deal(false);
        [b.qsmooth] = deal([]);
        [b.xsmooth] = deal([]);
    end %clearsmoothing()
    
    
    % ------------------------------------------
    %function tf = isregistered
    %    % May be coming soon...
    %    keyboard   
    %end
    
    
    % ------------------------------------------
    function b = smoothpose(b,sfun)
        %sfun = @(y)smooth(theta,y,width,'rloess');
        
        n = numel(b);
        for j = 1:n     % ==> upgrate to parfor
            if ~isempty(b(j).qraw) && ~all(isnan(b(j).qraw))
                qj = b(j).qraw;   %\_ raw data
                xj = b(j).xraw;   %/
                
                qj = qj.smooth(sfun,1);         %- Smooth q
                for c = 1:3                     %\ Smooth x
                    xj(:,c) = sfun(xj(:,c));    %/
                end
                
                b(j).qsmooth = qj;      % update 
                b(j).xsmooth = xj;      %
                b(j).smoothed = true;
            end
        end
        
    end %smoothpose()
               
    
    %=====================================================================
    % Getter Methods
    
    % ------------------------------------------
    function q = get.q(b)
        % Getter method for selecting the most processed state of q
        % - If a smoothed quaternion state exists, return that
        % - If only a raw quaternion state exists, return that
        % See also get.x()
        if ~isempty(b.qsmooth)
            q = b.qsmooth;
        else
            q = b.qraw;
        end
    end %get.q()
    
    % ------------------------------------------
    function x = get.x(b)
        % Getter method for selecting the most processed state of x
        % - If a smoothed position state exists, return that
        % - If only a raw position state exists, return that
        % See also get.q()
        if ~isempty(b.xsmooth)
            x = b.xsmooth;
        else
            x = b.xraw;
        end
    end %get.x()
    
    
    %=====================================================================
    % Setter Methods
    
    % ------------------------------------------
    function b = set.q(b,q)
        % Setter method for setting the state of q
        % q is a dependent property containing the rotation quaternion, so
        % when setting q, we are actually setting qraw and then clearing
        % qsmooth 
        b.qraw = q;
        b.qsmooth = [];
    end %set.q()
    
    % ------------------------------------------
    function b = set.x(b,x)
        % Setter method for setting the state of q
        % x is a dependent property containing the displacement vector so
        % when setting x, we are actually setting xraw and then clearing
        % xsmooth 
        b.xraw = x;
        b.xsmooth = [];
    end %set.x()
    
end %methods

end %classdef


% ------------------------------------------------------------------------
function b = struct2obj(s)
% Convert structure to bone object
%
% This function is primarily used for upgrading old structures to a
% current Bone class object
b = Bone();
b = repmat(b,size(s));

% These *will* exist:
[b.Tag] = deal(s.Tag);
[b.HiRes] = deal(s.HiRes);
[b.LoRes] = deal(s.LoRes);

% The pose states need more careful attention:
if isfield(s,'qraw')
    % This field will exist when some smoothing has been done
    [b.qraw] = deal(s.qraw);    
    [b.xraw] = deal(s.xraw);
    [b.qsmooth] = deal(s.q);
    [b.xsmooth] = deal(s.x);
else
    % Otherwise no smoothing has been done and only q & x exist
    [b.qraw] = deal(s.q);
    [b.xraw] = deal(s.x);
end

end %struct2obj()

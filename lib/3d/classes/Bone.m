classdef Bone
% BONE class.  
%
% Internally manages the raw/smoothed states of q & x
%
% Supersedes "Models" structure.

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
end

properties (Hidden=true, SetAccess=private) 
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
    %    keyboard
    %end
    
    
    % ------------------------------------------
    function b = smoothpose(b,sfun)
        %sfun = @(y)smooth(theta,y,width,'rloess');
        
        n = numel(b);
        parfor j = 1:n     % ==> upgrate to parfor
            if ~isempty(b(j).qraw) && ~all(isnan(b(j).qraw))
                qj = b(j).qraw.unwrap;   %\_ raw data
                xj = b(j).xraw;          %/
                
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

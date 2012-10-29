function varargout = ColourSpec(varargin)
%COLOURSPEC Return colour specifications for 8 standard Matlab colours.
%
% This function uses British spelling in its name primarily because this is
% the correct spelling for "colour", but also to avoid conflict with the
% Matlab help topic 'colorspec'.
%
% MAP = COLOURSPEC() returns an 8-by-2 cell array containing standard colour
% names in the first column, and their corresponding RGB colour triplets in
% the second column
%
% RGB = COLOURSPEC(COLOURNAMES) takes in N-by-1 cell array of colour names
% and returns an N-by-3 matrix of corresponding RGB colours.  COLOURNAMES
% can contain short names, ie, 'w', 'k', 'b', etc.
%
% COLOURNAMES = COLOURSPEC(RGB) is the reverse of the above, and takes in an
% N-by-3 matrix of RGB colours and returns their corresponding colour names
% in COLOURNAMES.  If any row in RGB is not a standard colour, it is
% returned as 'black'.
%
% RGB = COLORSPEC('TORGB',MIXED_CELL) converts a mixed cell of rgb colours
% and/or short or long names into RGB colour triplets
%
% COLOURNAMES = COLOURSPEC('TOLONGNAMES',MIXED_CELL) converts a mixed cell of
% rgb colours and/or short names and/or long names to long names.
%
% Only minmal error checking.
%
%   See also LINESPEC.

% Joshua Martin, 10-Feb-2012

% Colour definitions
map = {...
    'blue',    'b',  [0 0 1];
    'green',   'g',  [0 1 0];
    'red',     'r',  [1 0 0];
    'magenta', 'm',  [1 0 1];
    'cyan',    'c',  [0 1 1];
    'yellow',  'y',  [1 1 0];
    'black',   'k',  [0 0 0];
    'white',   'w',  [1 1 1]};

if nargin == 0
    varargout{1} = map;
    return
end

opt = 'auto';

if nargin == 1;
    inpt = varargin{1};
elseif nargin == 2;
    [opt, inpt] = varargin{:};
end


switch lower(opt)
    case 'auto'
        if iscell(inpt) && ischar(inpt{1})
            % COLORSPEC(COLORNAMES)
            % List of colour names provided - return RGB
            varargout{1} = mixed2rgb(inpt);
            
            
        elseif isnumeric(inpt) && (size(inpt,2) == 3)
            % COLORSPEC(RGB)
            % Matrix of RGB triplets provided, return colour names
            varargout{1} = rgb2names(inpt);
                        
        else
            error('Bad input argument');
        end
        
    case 'torgb'
        % COLORSPEC('ToRGB',MIXED_CELL)
        % MIXED_CELL is a cell array which could contain rgb triplets,
        % short names, or longnames
        
        varargout{1} = mixed2rgb(inpt);
        
    case 'tolongnames'
        % COLORSPEC('ToLongnames',MIXED_CELL)
        
        varargout{1} = mixed2longnames(inpt);
        
end %switch


% ------------------------------------------------------------------------
    function longnames = mixed2longnames(inpt)
        rgb = mixed2rgb(inpt);
        longnames = rgb2names(rgb);
        longnames = longnames(:)';
    end %mixed2longmanes()

%  ------------------------------------------------------------------------
    function rgb = mixed2rgb(mixed)
        if isnumeric(mixed)
            rgb = mixed;
            return
        elseif iscell(mixed)
            rgb = zeros(numel(mixed),3);    % default is black
            for j = 1:numel(mixed)
                cj = mixed{j};
                if isnumeric(cj)
                    % RGB triplets get injected directly
                    rgb(j,:) = cj(:)';
                else
                    % Strings get looked up
                    tf = strcmpi(map,cj);
                    if any(tf(:))
                        [row,~] = find(tf);
                        rgb(j,:) = map{row,end};
                    end
                end
            end
        else
            error('Bad input type')
        end
        
    end

% ------------------------------------------------------------------------
    function cnames = rgb2names(RGB)
            cnames = cell( size(inpt,1), 1 );
            rgbList = cell2mat(map(:,3));
            
            % Populate list with names:
            for j = 1:size(RGB,1)
                idx = all( rgbList == repmat( RGB(j,:), size(rgbList,1),1 ) ,2);
                if any( idx )
                    cnames{j} = map{idx,1};     % Assign colour name
                else
                    cnames{j} = 'black';        % Default
                end
            end
    end %rgb2names()

end %ColorSpec()

function kpf2d(hfig,event)

% Key listings:
%   leftarrow    Cycle phases down
%   rightarrow   Cycle phases up
%   uparrow      Cycle slices up
%   downarrow    Cycle slices down
%   
%   shift+leftarrow   Pan left
%   shift+rightarrow  Pan right
%   shift+uparrow     Pan up
%   shift+downarrow    Pan down
%               '+'    Zoom in
%               '-'    Zoom out

% First we set the appdata 'CurrentKey' so we can determine at any point
% within the program what key is currently depressed:
setappdata(hfig,'CurrentKey',event.Key)

% Now manage the key presses:
handles = guidata(hfig);

if isempty(handles.DICOM) || isempty(handles.DICOM.X)
    return
end

ha = handles.axes1;         % Shorthand
s = current('slice',ha);    
p = current('phase',ha);

X = handles.DICOM.X;        % Shorthand
ns = size(X,3);
np = size(X,4);

switch event.Key
    
    case {'uparrow','downarrow','leftarrow','rightarrow'}
        
        % Cycle through the image stack with arrow keys:
        ds = 0;
        dp = 0;
        if isempty(event.Modifier)
            switch event.Key
                case 'uparrow'
                    ds = +1;
                case 'downarrow'
                    ds = -1;
                case 'rightarrow'
                    dp = +1;
                case 'leftarrow'
                    dp = -1;
                otherwise
                    return
            end
            
            s = s+ds;                   % Incrment slice
            s = max( 1, min(ns,s) );    % apply bounds
            
            p = p+dp;                   % Increment phase
            p = max( 1, min(np,p) );    % apply bounds
            
            updateSlice(handles,s,p)
            
            
        % Pan the image around with the shift key
        elseif any(strcmpi(event.Modifier,'shift'))
            %keyboard
            x = get(ha,'XLim');
            y = get(ha,'YLim');
            d = round(0.08*diff(x));        % 20% in each step
            switch event.Key
                case 'uparrow'
                    y = panlims(-d,y,size(X,1));
                case 'downarrow'
                    y = panlims(d,y,size(X,1));
                case 'rightarrow'
                    x = panlims(d,x,size(X,2));
                case 'leftarrow'
                    x = panlims(-d,x,size(X,2));
            end
            %keyboard
            set(ha,'XLim',x)
            set(ha,'YLim',y)
            
        end
        
        
    case {'equal','hyphen','add','subtract'}
        % User pressed + or - for zooming:
        s = 1.5;
        switch event.Character
            case '+' %{'=','+'}
                s = s; % zoom in
            case '-' %{'-','_'}
                % zoom out
                s = 1/s;
            otherwise
                return
        end
        
        zoom(ha,s);
        
        
    otherwise
        %disp(event)
        %keyboard

end


% ------------------------------------------------------------------------
function lims = panlims(d,lims,imdim)
%      d: 1-by-1 defining how much we are to shift by
%   lims: 1-by-2 of the current axis limits to be modified
% imdims: 1-by-1 dimension of the image in the requested orientation

% Normal axis limits for an image of 100x100 are [0.5 100.5], so 1/2 pix
% each side of the numerical dimensions

% Add a buffer around the image that we can pan to:
buff = 30;

% limit motion to image bounds:
dmax = imdim - (lims(2)-0.5) + buff;  % Max allowable positive shift
dmin = 1 - (lims(1)+0.5) - buff;      % Max allowable negative shift

% Apply limits to the diplacement:
d = min([dmax d]);
d = max([dmin d]);

% Adjust axis limits
lims = lims+d;

% Conditinal for zooming out past fit - zoom with image at centre:


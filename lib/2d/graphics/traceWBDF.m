function traceWBDF(hf,~)
%TRACEWBDF WindowButtonDownFcn callback which instigates a tracing
%
% See also TRACEWBMF, TRACEWBUF (nested functions)

% This function will capture all mouse clicks, but we're only actioning on
% normal clicks: 
if ~strcmp('normal', get(hf,'SelectionType'))
    return
end

% Graphics pad on windows instigates a trace using a 'normal' selection,
% but then has the possibility of instigating another 'normal' selection
% using the button on the pen (or conversely if started with pen button,
% then second click can come from placing the pen on the pad).  We avoid
% this by checking if a trace is already in progress
if ~isempty( get(hf,'WindowButtonUpFcn') )
    % Tracing still underway; bail out
    return
end

% Get handles:
handles = guidata(hf);
ax = handles.axes1;

% Discretisation control:
%   Specify min & max spacing for points on a trace:
DTGT = 1;%mean(handles.Images.s); % Target discretisation
DMIN = DTGT/2;
DMAX = DTGT*4;

% Joining tolerance:
%   If tracing within this distance (in pixels) to an existing roi, it 
%   will join with that roi:
JOINTOL = 3;

% Get the initial point
ptmat = get(ax,'CurrentPoint');
xs = ptmat(1,1);
ys = ptmat(1,2);

% Create a new working trace object:
tobj = line(xs,ys,'tag','tmpTrace','Color','r');

% See if we should be be joining to a parent
rparent = [];
get_parent_roi()


% Add the callbacks for tracing:
oldwbmf = get(hf,'WindowButtonMotionFcn');
set(hf,'WindowButtonMotionFcn',@traceWBMF)
set(hf,'WindowButtonUpFcn', @traceWBUF);

% Now the window button motion function (TRACEWBMF) takes control

    % --------------------------------------------------------------------
    function traceWBMF(~,~)
        %TRACEWBMF WindowButtonMotionFcn callback which draws a trace during creation
        %
        % See also TRACEWBDF, TRACEWBUF
        
        tmpx = get(tobj,'xdata');
        tmpy = get(tobj,'ydata');
        a=get(gca,'CurrentPoint');
        p = a(1,1:2);
        
        % Discretisation control:
        %   Only current point to trace if it is further than DMIN from the
        %   last trace point:
        last = [tmpx(end), tmpy(end)];
        dc = p - last;
        d = sqrt(sum(dc.*dc));
        if (d >= DMIN) && (d <= DMAX)
            tmpx(end+1) = p(1);
            tmpy(end+1) = p(2);
            set(tobj,'xdata',tmpx,'ydata',tmpy);
        elseif d > DMAX
            % If it is further than DMAX, add in some points:
            n = ceil(d/DMAX);
            xv = linspace(tmpx(end),p(1),n+1);
            yv = linspace(tmpy(end),p(2),n+1);
            tmpx(end+1:end+n) = xv(2:end);
            tmpy(end+1:end+n) = yv(2:end);
            set(tobj,'xdata',tmpx,'ydata',tmpy);
        end
        
        % Update the annotations
        %   (This doesn't happen automatically because we removed the old wbmf)
        infoWBMF(hf,[])
        
    end

    % --------------------------------------------------------------------
    function traceWBUF(~,~)
        %TRACEWBUF WindowButtonUpFcn callback which finishes a tracing
        %
        % hf is the handle to the current figure
        %
        % See also TRACEWBDF, TRACEWBMF
        
        
        % Get Currently displayed slice:
        s = current('slice',handles.axes1);
        p = current('phase',handles.axes1);
        
        % Get the full path of that Dicom file:
        dcm = [handles.Images.pth, handles.Images.files{s,p}];
        
        % Get its header data:
        PS  = handles.Images.info(s,p).PixelSpacing;
        IPP = handles.Images.info(s,p).ImagePositionPatient;
        IOP = handles.Images.info(s,p).ImageOrientationPatient;
        
        % Clear the drawing callbacks:
        set(hf,'WindowButtonMotionFcn',[])
        set(hf,'WindowButtonUpFcn',    [])
        
        % Get line data
        x = get(tobj,'xdata')';
        y = get(tobj,'ydata')';
        delete(tobj)                 % Delete object
        
        % Curves with less than 5 points are ignored:
        if numel(x) < 5
            % do nothing
            
        else
            % Switch if we're creating a new one or joining to exisitng:
            if isempty(rparent)
                                
                % Get the default colour:
                kids = get(handles.MI_DefROIColour,'Children');
                clr = get(kids(strcmpi('on',get(kids,'Checked'))),'Label');
                
                % Now create new roi
                r = roi(x,y,s,p,dcm,PS,IPP,IOP,clr);
                r = r.equispace(DTGT);
                r = r.close;
                handles.traces(end+1) = r;
            else
                % Join to roi recorded as parent
                ind = strcmp({handles.traces.Tag},rparent.Tag);
                r = rparent.join([x(:) y(:)]);
                r = r.equispace(DTGT);
                handles.traces(ind) = r;
            end
            
            % Update handles - this must come before updateSlices() because
            % that function triggers other functions which are dependent on
            % handles being up to date.
            guidata(hf,handles)
            
        end
        
        % Refresh traces from handles for current slice/phase
        %   We do this here so that if user began editing an roi, but the
        %   edit had less than 5 points and didn't get drawn, the parent
        %   roi would still be bold.  But calling updateSlice will refresh
        %   all the drawing.
        updateSlice(handles);
        
        % Restore old WBMF
        set(hf,'WindowButtonMotionFcn',oldwbmf)
           
    end %traceWBUF()
    
    
    % --------------------------------------------------------------------
    function get_parent_roi
        if isempty(handles.traces)
            return
        end
        % Get Currently displayed slice:
        s = current('slice',handles.axes1);
        p = current('phase',handles.axes1);
        
        % Select rois on this slice / phase
        inds = ([handles.traces.Slice] == s) & ([handles.traces.Phase] == p);
        rset = handles.traces(inds);
        
        % Find out how far away they are:
        n = numel(rset);
        d = NaN(1,n);
        for j = 1:n
            d(j) = min( hypot( (xs - rset(j).x), (ys - rset(j).y) ) );
        end
        
        % Choose the closest roi & get its distance
        [d,j] = min(d);
        r = rset(j);
        
        % If that roi is closer than the join tolerance, we're going to
        % modify it:
        if d <= JOINTOL
            rparent = r;    % Store for modifying on exit
            hr = findobj(handles.axes1,'Tag',rparent.Tag);
            set(hr,'LineWidth',1.5)
        end
        
    end %get_trace_parent()
    


end
function roi_user_interaction_test

% Get the gui and make sure it has an image loaded:
hf = prep_gui_2D;

% These can be a bit unstable on OSX because sometimes Matlab doesn't
% refresh its graphics properly when maximising/restoring windows...  In
% that case, restart Matlab and try again.
test_drawing_to_scale(hf);
test_joining(hf)


%error('More tests to follow')
% #Importing
% #Exporting
% #Navigator

%hf = prep_gui_3D;
% #3D view
% #Loading into 3D workspace
%


end


% =========================================================================
% TEST FUNCTIONS:


% ------------------------------------------------------------------------
function test_drawing_to_scale(hf)

% Use java to draw the roi:
[xy,L,A] = new_trace_path('AreaTest');
jdraw(hf,xy);

% Examine trace:
handles = guidata(hf); 
% figure, plot(handles.traces.x,handles.traces.y,'r')
% hold on
% htrue = plot(xy(:,1),xy(:,2),'b');
% copyobj(htrue,handles.axes1)
r = handles.traces(end);

% Check that the area and length are correct:
assertAlmostEqual(r.Area.pixels,   A, 1/100)    % 1% for drawing tolerance
assertAlmostEqual(r.Length.pixels, L, 1/100)    % 1% for drawing tolerance

end %test_drawing_to_scale()



% ------------------------------------------------------------------------
function test_joining(hf)
% This function draws two rois and checks that traceWBDF>traceWBUF makes
% correct calls to join two rois when required.

% TEST 1 - Join when starting near an exisiting roi:

[xy1,xy2,xyResult] = new_trace_path('JoinTest');

jdraw(hf,xy1);

jdraw(hf,xy2);


% TEST 2 - Do not join when not close enough to exisiting roi:

end









% =========================================================================
% =========================================================================
% =========================================================================
% Helper functions


% ------------------------------------------------------------------------
function jdraw(hf,xy)
% Use a java robot to draw something, then, hand back out the trace object
handles = guidata(hf);
nt = numel(handles.traces);
ha = handles.axes1;

% Here's something really annoying
%  set(0,'HideUndocumentd','off') % show properties
wtf = get(ha,'WarpToFill');
set(ha,'WarpToFill','on')

% We'll be more accurate if we make the plot bigger:
fpsn = get(hf,'Position');
jf = getjframe(hf);
jf.setMaximized(true);
drawnow, pause(0.05);

% Turn the tracing tool on:
set(handles.TraceTool,'State','on')
Segmentation('TraceTool_ClickedCallback',handles.TraceTool,[],handles)

% Convert axes coords to screen coords:
xys = axisloc2screenloc(ha,xy);

% Convert to java screen space for the robot:
xyj = matlabpsn2javapsn(xys);

% Instigate a robot:
%   Java positions are (0,0) in the top left.  
mouse = java.awt.Robot;
btn1 = java.awt.event.InputEvent.BUTTON1_MASK;
%btn2 = java.awt.event.InputEvent.BUTTON2_MASK;

% Now draw:
figure(handles.figure1), drawnow    % necessary
pause(0.05)
dt = 0.1;
for j = 1:size(xyj,1) 
    mouse.mouseMove(xyj(j,1),xyj(j,2));
    drawnow, pause(dt) % so our drawing function can try to keep up
    if j == 1
        mouse.mousePress(btn1)
        drawnow, pause(dt)
    end
end 
mouse.mouseRelease(btn1);
drawnow

% Now wait for handles to update:
t = tic;
while numel(handles.traces) == nt && toc(t) < 0.5
    handles = guidata(hf);
    pause(0.05)
end
% ...and we're done, except for this silly thing:
set(ha,'WarpToFill',wtf)

% Oh, and this
jf.setMaximized(false)
set(hf,'Position',fpsn)

end %jdraw()


% ------------------------------------------------------------------------
function h = new_gui_handle()
h = findall(0,'Name','Segmentation');
try close(h), end %#ok<TRYNC>
h = Segmentation; 
drawnow
end %get_gui_handle()


% ------------------------------------------------------------------------
function varargout = new_trace_path(tag)

switch tag
    case 'AreaTest'
        % [xy,L,A] = new_trace_path('AreaTest')
        x = [0 2 2 3 3 5 5 7 7 5 5 4 4 1 1 0 0]*20 + 100;
        y = [2 2 0 0 4 4 3 3 6 6 7 7 6 6 3 3 2]*20 + 100;
        xy = [x(:) y(:)];
        L = sum( sqrt( diff(x).^2 + diff(y).^2) );
        A = polyarea(x,y);
        varargout = {xy,L,A};
        
    case 'JoinTest'
        % [xy1,xy2,xyJoined] = new_trace_path('AreaTest')
        
        th = linspace(0,2*pi,100)';
        xy1 = [cos(th) sin(th)];
        
        th = linspace(pi-pi/6,pi+pi/10,20)';
        xy2 = [cos(th)*2 + 1.5,  sin(th)*2];
        
        % Re-scale:
        s = 80;
        t = [150 150];
        xy1 = xy1*s + repmat(t, size(xy1,1),1);
        xy2 = xy2*s + repmat(t, size(xy2,1),1);
        
        % Use the roi method to do the join:
        r1 = roi(xy1(:,1),xy1(:,2));
        r2 = r1.join(xy2);
        xyResult = [r2.x, r2.y];
        
        % Hand out variables:
        varargout = {xy1,xy2,xyResult};
        
end


end %new_trace_path()


% ------------------------------------------------------------------------
function h = prep_gui_2D

h =  new_gui_handle();

handles = guidata(h);

% Load an image set
disp('Loading data ...')
%saved_data = load('simple_dicom_set.mat');
%handles.DICOM = saved_data.DICOM;
%disp('Complete!')

hObject = h; % For convenience

%
% Specify an MATLAB example image:
pathname = [toolboxdir('images') filesep 'imdemos' filesep];
filename = 'knee1.dcm';

% The following code is ripped straight out of
% Segmentation>MI_LoadDicom_Callback()
handles.DICOM.pth = pathname;
if ~iscell(filename), filename = {filename}; end
try
    [X,z,s,dinfo] = readDicoms( pathname, filename );
catch ME
    if strcmp(ME.identifier,'MATLAB:waitbar:InvalidInputs')
        disp('User cancelled')
        return
    else
        rethrow(ME)
    end
end
handles.DICOM.files = filename(:);  % Slices populate a column
handles.DICOM.info  = dinfo;
handles.DICOM.X   = X;
handles.DICOM.z   = z;
handles.DICOM.s   = s;
handles.DICOM.CLim = imlimits(handles.DICOM.X,0.9999);

% Update user path:
handles.DICOM.pth = pathname;
handles.userPath = pathname;
%}

% Update guidata
guidata(hObject,handles)

% Now that the guidata is updated, configure the view:
configureView2(handles);
set(handles.axes1,'CLim',handles.DICOM.CLim)

   
% Now we have the handle, and a DICOM image is loaded.
    % ------------------------------
    function clim = imlimits(X,frac)
        % See the code of STRETCHLIM.  It returns the limits in percentage of the
        % data range, so we need to multiply this by the number of bins to return
        % the actual intensity value
        if isa(X,'uint8')
            nbins = 256;
        else
            nbins = 65536;
        end
        lowhigh = stretchlim(X(:),[0 frac]); % discard 0.01% of outliers;
        clim = lowhigh*nbins;
        clim = clim(:)';
    end

end %prep_gui_2D()





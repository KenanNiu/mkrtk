classdef FigLocker < handle
% FigLocker Class for bocking figures
%
%

properties
    fig
end

% State property
properties (SetAccess=protected)
    islocked = false;
end

% Target/host properties
properties (Hidden=true)
    Resize
    figjFrame
    glassJWindow
    controlJWindow
    psnCbks
    listeners
    barData 
    labelData
end

% Private properties
properties (Hidden=true,SetAccess=protected,GetAccess=protected)
    jProgressBar
    jLabel
    jButtons
end

% Normal methods
methods
    
    % ----------------------------------------------
    function h = FigLocker( hf )
        h.listeners = handle([]);   % Ideally this would be set as a default
                                    % in the properties list above, but
                                    % that causes problems with object 
                                    % disposal.  So we set it here.
        switch nargin
            case 0
                % Empty locker
                return
                
            case 1
                % Get or instantiate locker:
                if isFigLocked(hf)
                    h = getlocker( hf );
                else
                    h.fig = hf;
                    h.Resize = get(hf,'Resize');
                    % The following listener ensures that the locker object
                    % gets destroyed correctly when a figure gets closed.
                    h.listeners(end+1,1) = ...
                        addlistener(hf,'ObjectBeingDestroyed',@destroyLocker);
                end
            otherwise
                error('Incorrect number of input arguments')
        end %switch   
    end %FigLocker()
    
    % ----------------------------------------------
    function addbuttons(h, varargin)
        assert(numel(varargin)>=2,...
            'A string and callback must be specified for each button')
        txt = varargin(1:2:end-1);
        cbk = varargin(2:2:end);
        [width,height] = calcMinButtonDims(txt);
        gap = 10;
        
        % If buttons already exist, we might want to delete them.  How?
        % removeComponent?
        
        % Get or instantiate a control window
        configureControlWindow(h);
        
        % Configure buttons
        wdim = h.controlJWindow.getSize;
        W = wdim.getWidth;
        H = wdim.getHeight;
        
        n = numel(txt);
        for j = n:-1:1
            jgaps = j - (n-1)/2 - 1;
            jwidths = j - n/2 -1;
            jb(j) = javax.swing.JButton(txt{j});
            h.controlJWindow.getLayeredPane.add(jb(j));
            jb(j).setSize(width,height);
            x = W/2+jgaps*gap+jwidths*width;
            y = H/2+0.8*height;
            jb(j).setLocation(x,y)
            jbtn = handle(jb(j),'CallbackProperties');
            set(jbtn,'ActionPerformedCallback',cbk{j})
        end
        h.jButtons = jb;
    end %addbuttons()
    
    % ----------------------------------------------
    function jProgbar = addprogressbar( h , indeterminate, value)
        if ~h.islocked
            error('This figure is not locked')
        end       
        
        % Set up properties:
        if nargin >= 2
            s.Indeterminate = indeterminate;
            if nargin >= 3
                s.Value = value;
            end
        end
        
        % Get or instantiate a control window
        configureControlWindow(h);
        
        % Make dimensions for bar:
        rect = h.controlJWindow.bounds;
        ww = rect.getWidth;
        wh = rect.getHeight;
        bw = h.TOOLWIDTH();
        bh = 20;
        bx = floor( ww/2 - bw/2 );
        by = wh/2-bh/4;
        bounds = java.awt.Rectangle(bx,by,bw,bh);
                
        % Instigate a new progress bar:
        s.Bounds = bounds;
        jProgbar = newProgbarWithProps(s,h.controlJWindow);
                
        % Store its handle
        h.jProgressBar = jProgbar;
        jProgbar.setBounds(bounds)
        
    end %addprogressbar()
    
    % ----------------------------------------------
    function settext( h, string )
        if ~h.islocked
            error('This figure is not locked')
        end
        
        % Get or instantiate a control window
        configureControlWindow(h);
        
        % Create new label if needed:
        if isempty(h.jLabel)
            h.jLabel = newLabel(h.controlJWindow);
        end
        
        % Set/update text:
        h.updateLabelString(string)
        
    end %settext()
    
    % ----------------------------------------------
    function setprogress( h, prog )
        assert(isnumeric(prog),'Progress should be numeric: 0->100, Inf, or NaN')
        if ~h.islocked
            error('This figure is not locked')
        end
        
        % Configure indeterminate / value
        if isfinite(prog)
            INDET = false;
        else
            INDET = true;
            prog = 0;
        end
        
        % Check if progress bar exists:
        if isempty(h.jProgressBar)
            h.addprogressbar(INDET,prog);
        else
            h.jProgressBar.setIndeterminate(INDET);
            if ~INDET
                h.jProgressBar.setValue(prog);
            end
        end
        
    end %setprogress()
    
    % ----------------------------------------------
    function lock( h, hf )
        % lockfig(h)        if h.fig already specified
        % lockfig(h,hf)     
        if nargin > 1 || isempty( h.fig )
            h.fig = hf;
        end
        if h.islocked
            return
        end         
        
        % Flag:
        h.islocked = true;
        
        figure(h.fig)   % ensure it's at the front
        drawnow         % just to ensure that hf is showing
        
        % Store locker's handle in the figure's appdata:
        setappdata(h.fig,TAG,h)
        
        % Get the figure's jframe
        %   Actually jFrame has class:
        %       com.mathworks.hg.peer.FigureFrameProxy$FigureFrame
        h.figjFrame = getJFrame(h.fig);
        
        % Disable figure
        h.figjFrame.setEnabled(false);
        
        % Create locking window & callbacks
        [h.glassJWindow,h.psnCbks{1}] = createGlassJWindow( h.figjFrame );
        
        % Callbacks & listeners
        callbackHelper('add',h,@(a,b)psnUpdateFcn(a,b,h))
        h.listeners(end+1,1) = addlistener(h.fig,'Visible','PostSet',@visHelper);
        
        % Prevent resizing
        h.Resize = get(h.fig,'Resize');     % Get state
        %set(h.fig,'Resize','off')           % Set off
             
    end %lock()
    
    % ----------------------------------------------
    function unlock( h )
        % Unlock target figure
        
        % For step numbers, see lock()
        
        % (4) Enable figure
        h.figjFrame.setEnabled(true);
           
        % (3) Restore resizing
        set(h.fig,'Resize',h.Resize)
        
        % (2) Remove callbacks & listeners
        callbackHelper('remove',h)
        h.listeners( ishandle(h.listeners) ).delete;    % Delete valid listeners
        
        % (1) Remove blocking jFrame if locked
        %   Force everything else (ie, not conditional)
        %   in case the glassJWindow somehow got closed and 
        %   we didn't catch it.  But we can only dispose
        %   of the glassJWindow if we have a handle for it, so
        %   we check with islocked()
        if isFigLocked( h.fig )            
            h.glassJWindow.dispose();
            if ~isempty(h.controlJWindow)
                h.controlJWindow.dispose();
            end
        end
        
        % (0) Remove object from appdata
        if isappdata(h.fig,TAG)
            rmappdata(h.fig,TAG)
        end
        
        % Clean up:
        %   Leave the following in place in case we want
        %   to lock again using this object:
        %       h.figjFrame 
        %       h.Resize
        h.controlJWindow = [];
        h.glassJWindow = [];
        h.jLabel  = [];
        h.jProgressBar = [];
        h.listeners(:) = [];
        
        % Flag: 
        h.islocked = false;
        
    end %unlock()
    
    % ----------------------------------------------
    function delete( h )
        % Called when the object gets cleared from memory.
        % If the figure is still locked, unlock it, then the instance will
        % get disposed of when the handle goes out of scope
        if isvalid(h)
            h.unlock
        end
    end %delete()
    
            
end %methods

% Static methods
methods (Static)
    
    % ----------------------------------------------
    function h = GetLocker( hf )
        h = getlocker( hf );
    end %GetLocker()
    
    % ----------------------------------------------
    function tf = IsLocked( hf )
        % Set to see if this figure is locked
        h = getlocker( hf );
        tf = ~isempty( h ) && h.islocked;
    end %IsLocked()
    
    % ----------------------------------------------
    function h = Lock( hf )
        h = FigLocker(hf);
        if ~h.islocked
            h.lock;
        end
    end %Lock()
    
    % ----------------------------------------------
    function Unlock( hf )
        h = getlocker( hf );
        if ~isempty(h)
            h.unlock % h.delete
        end
    end %UnLock()
    
end %methods (Static)

% Private methods
methods (Access=private)
    % ----------------------------------------------
    function dohide(h)
        % We need this function to handle the blocking window when the host
        % figure gets hidden.
        % Similar to unlock() but for use when a locked figure gets
        % hidden
        
        % If a progbar exists, store its state data:
        if ~isempty(h.jProgressBar)
            h.barData = setgetProgbarData('get',h.jProgressBar);
            h.jProgressBar = []; % drop handle
        end
        % If text exists, store its data:
        if ~isempty(h.jLabel)
            h.labelData = setgetLabelData('get',h.jLabel);
        end
        h.glassJWindow.dispose;
        h.glassJWindow = [];
        if ~isempty(h.controlJWindow)
            h.controlJWindow.dispose;
            h.controlJWindow = [];
        end
        %h.unlock;     % Unlock the figure, it will get locked again by the doshow
    end
    
    % ----------------------------------------------
    function doshow(h,varargin)
        % Similar to lock() but for use in a listener on 'Visible' property
        % when a hidden locked figure gets made visible.  It locks the
        % figure (because this function will only be called when a locked
        % figure is becoming visible) and restores the children, if any.
        %
        % varargin exists to keep the invoking timer from complaining.
        h.islocked = false;
        h.lock;
        % Create control window if necessary
        if ~isempty(h.barData) || ~isempty(h.labelData)
            configureControlWindow(h)
        end
        % Re-create progbar if one existed:
        if ~isempty(h.barData)
            % Make a new bar
            h.jProgressBar = newProgbarWithProps(h.barData,h.controlJWindow);
            % Clear data:
            h.barData = [];
        end        
        % Re-create text item, if one existed
        if ~isempty(h.labelData)
            % Make new text:
            h.jLabel = newLabel(h.controlJWindow);
            h.updateLabelString(h.labelData.Text.toString)
            % Clear data:
            h.labelData = [];
        end
        
        % Restore callback(s)
        %callbackHelper('add',h,h.psnCbk)
    end
    
    % ----------------------------------------------
    function configureControlWindow(h)
        if isempty(h.controlJWindow)
            [h.controlJWindow,h.psnCbks{end+1}] = createControlJWindow(h.figjFrame);
        end
    end %configureControlWindow()
    
    % ----------------------------------------------
    function w = TOOLWIDTH(h)
        % Convenient central place for calculating/defining the (maximum) width of
        % the text/progress bar:
        rect = h.glassJWindow.bounds;
        ww = rect.getWidth;
        w = floor(ww*5/8);
    end %WIDTH()
    
    % ----------------------------------------------
    function updateLabelString(h,str)
        
        % Cells   -> format as is
        % Strings -> wrap to width
        if iscell(str)
            C = str(:)';
            C(2,:) = {sprintf('\n')};
            C = C(:); C(end) = [];
            str = [C{:}];
        else
            % Need to implement text wrapping for strings here...
        end
        
        % Current dimensions for text:
        siz = h.jLabel.getSize;
        w_old = siz.width;
        h_old = siz.height;
        
        % Update text & get its new dimensions
        r0 = java.awt.Rectangle(0,0,0,0);
        h.jLabel.setText(str);  % This causes it to expand
        h.jLabel.setBounds(r0); % So shrink it quickly before it gets drawn, just while we're still working on it...
        
        % New preferred dims:
        siz = h.jLabel.getPreferredSize;
        w_new = siz.width + 8;
        h_new = siz.height + 3;
        
        w_max = h.TOOLWIDTH();
        
        % Make dimensions for text:
        rect = h.controlJWindow.bounds;
        ww = rect.getWidth;
        wh = rect.getHeight;
        th = 20;
        
        if w_new < w_old
            tw = w_old;
        elseif w_new > w_max
            tw = w_max;
        else
            tw = w_new;
        end
        
        
        tx = floor( ww/2 - tw/2 );
        ty = floor( wh/2 - th-10 );
        
        rect = java.awt.Rectangle(tx,ty,tw,th);
        
        h.jLabel.setBounds(rect);
        
        
    end %updateLabelString()
    
end %methods (Private)

end %classdef()



% ======================== Helper Functions ============================ %

% ------------------------------------------------------------------------
function [w,h] = calcMinButtonDims(textcell)
% Calculate the dimensions required which will be appropriate for all the
% strings in the cell.
jL = javax.swing.JButton;
for j = numel(textcell):-1:1
    jL.setText(textcell{j});
    dims = jL.getPreferredSize();
    w(j) = dims.getWidth;
    h(j) = dims.getHeight;
end
w = max(w);
h = max(h);
end %calcButtonDims()

% ------------------------------------------------------------------------
function jtxt = newLabel(jWindow)
% S must be at a minimum:
%   s.Bounds    => java.awt.Rectangle
%   s.Text      => string

% Create jlabel
jtxt = javax.swing.JLabel;

% Set properties
jtxt.setHorizontalAlignment(0) %centred
jtxt.setVerticalAlignment(0)

% Add to Layered pane:
jWindow.getLayeredPane.add( jtxt );
jtxt.show();

end %newLabelWithProps()


% ------------------------------------------------------------------------
function jpb = newProgbarWithProps(s,jWindow)
assert(isfield(s,'Bounds'))
% Defaults:
if ~isfield(s,'Value')
    s.Value = 0;
end
if ~isfield(s,'String')
    s.String = java.lang.String([num2str(s.Value),'%']);
end
if ~isfield(s,'Indeterminate')
    s.Indeterminate = true;
end
% Create a progress bar:
jpb = javax.swing.JProgressBar;

% Set Properties
if ~s.Indeterminate
    jpb.setValue(s.Value)
end
jpb.setIndeterminate(s.Indeterminate)
jpb.setBounds(s.Bounds)
jpb.setString(s.String)
jpb.setOpaque(true)
%jpb.setPreferredSize(jProgbar.getSize)

% Add it to the content pane:
jWindow.getLayeredPane.add( jpb );
jpb.show();
end %newProgbarWithProps()

% ------------------------------------------------------------------------
function callbackHelper(opt,h,psnCbk)
jFrameCbkProps = handle(h.figjFrame,'CallbackProperties');
switch opt
    case {'add','set'}
        % Set some callbacks
        set(jFrameCbkProps,'ComponentMovedCallback',psnCbk)
    case 'remove'
        % Remove the callbacks
        set(jFrameCbkProps,'ComponentMovedCallback','')
    otherwise
        error('bad option')
end
end %callbackHelper()

% ------------------------------------------------------------------------
function destroyLocker( hf, ~ )
% Clean up when the figure hf is being destroyed.
% Because we have a circular reference (hf stores h, and h stores hf), we
% need to unlink the two handles before destruction, otherwise matlab will
% hang on to an instance of the class.  Not terrible, but annoying when you
% try to do "clear classes" and it won't.  Do we do the following when the
% figure's ObjectBeingDestroyed event gets triggered:
hf = double(hf);
h = getlocker(hf);
rmappdata(hf,TAG)
delete(h)
clear h
end %destroyLocker()

% ------------------------------------------------------------------------
function tf = isFigLocked( hf )
h = getlocker(hf); 
if isempty(h) || isempty(h.glassJWindow) || ~h.glassJWindow.isDisplayable 
    tf = false;
else
    tf = true;
end
end %islocked()

% ------------------------------------------------------------------------
function h = getlocker( hf )
% Convenience wrapper for getappdata call
h = getappdata(hf,TAG);
end %getlocker()

% ------------------------------------------------------------------------
function t = TAG()
% Convenient place for defining the tag by which we store the locker in the
% host figure's appdata
% Used by:
%   getappdata(hfig,TAG)
%   setappdata(hfig,TAG)
    t = 'FigLocker_Object'; 
end %TAG()

% ------------------------------------------------------------------------
function [jWindow,psnfcn] = createControlJWindow(jFrame)
% Create a jWindow which is used to hold any controls and progress bar
% By the time we reach here, a glassJWindow will already exist
%
% See also createGlassJWindow

% Create the jWindow:
jWindow = javax.swing.JWindow(jFrame);
jWindow.pack();

psnfcn = @psnCtrlWindow;
psnCtrlWindow();

jWindow.show();

drawnow

    function psnCtrlWindow(varargin)
        if jFrame.isDisplayable
            jFrame.getGlassPane().setVisible(true); % Must be visible to get it's location
            fpsn = jFrame.getGlassPane().getLocationOnScreen;
            fsiz = jFrame.getGlassPane().getSize;
            wsiz = [fsiz.getWidth*2/3 120];
            x = fpsn.getX + fsiz.getWidth/2-wsiz(1)/2;
            y = fpsn.getY+fsiz.getHeight/2-wsiz(2)/2;
            wpsn = [x, y];
            
            wsiz = java.awt.Dimension(wsiz(1),wsiz(2));
            wpsn = java.awt.Point(wpsn(1),wpsn(2));
            
            jWindow.setLocation( wpsn );
            jWindow.setSize( wsiz );
            jFrame.getGlassPane().setVisible(false);
        end
    end

end

% ------------------------------------------------------------------------
function [jWindow,psnfcn] = createGlassJWindow(jFrame)
% hf - handle to the figure
% jFrame - the output of getJFrame(hf)
jFrame.show();              % ensure it's showing

% Create a transparent JWindow:
jWindow = javax.swing.JWindow(jFrame);
jWindow.pack();
t = tic;
while ~jWindow.isDisplayable() ...  % This ensures that the
        && toc(t) < 1.0             % jWindow is fully created
    disp('pausing')
    pause(0.02);                    % before moving on
end

% Set colour & transparency
%jWindow.setBackground( java.awt.Color(0.8,0.8,0.8,alpha) ); 
%".setBackground()" with transparency prevents correct refreshing of the bakcground!!
jWindow.setBackground( java.awt.Color(0.8,0.8,0.8) );
com.sun.awt.AWTUtilities.setWindowOpacity(jWindow, 0.75);
jWindow.show();

% Set the position/size of the jWindow & return function handle
psnfcn = @positionOver;
positionOver();

% Now we need to ensure that the window is displayed before exiting
% Seems like a simple drawnow call is what works.
drawnow


    % -------------------------------
    function positionOver(varargin)
        % There are occasional cases that slip through where this function
        % can be called after the jWindow has been destroyed.
        try
            jFrame.getGlassPane().setVisible(true); % Must be visible to get it's location
            jWindow.setLocation( ...                        % Set position of panel...
                jFrame.getGlassPane().getLocationOnScreen );% To the position of the glass panel
            jWindow.setSize(jFrame.getGlassPane().getSize); % Set size of panel to size of glass panel
            jFrame.getGlassPane().setVisible(false);    % Can turn this off again
        catch  %#ok<CTCH>
            % jFrame has probably been destroyed
        end
    end

end %createBlockingJframe()


% % ------------------------------------------------------------------------
% function psn = mpsn2jpsn(psn)
% ss = get(0,'ScreenSize');
% psn(1) = psn(1) - 1;                    % x-position
% psn(2) = ss(4) - (psn(2)-1) - psn(4);   % y-position
% end %mpsn2jpsn()

% ------------------------------------------------------------------------
function visHelper(~,eventdata)
% We could hide/show the blocking glassJWindow, but that causes problems with
% size & positioning of children when it gets re-created.  Better to just
% dispose of it, and re-create it when the figure becomes visible again.
% That means also re-creating any child components that have been added
% (progress bar, text)
hf = eventdata.Affectedobject;
h = getlocker(hf);
switch eventdata.NewValue
    case 'on' 
        start(timer('ExecutionMode','singleShot',...
            'StartDelay',0.05,'TasksToExecute',1,...
            'TimerFcn',@h.doshow,...
            'StopFcn',@(t,e)delete(t)))
    case 'off'
        h.dohide;
end        

end %visHelper()


% ------------------------------------------------------------------------
function varargout = setgetProgbarData(opt,varargin)
% Set/get selected progress bar properties
%   progbarData('set',jProgressBar,propstruct)
%   s = progbarData('get',jProgressBar)
switch opt
    case 'set'
        [jpb,s] = varargin{:};
        jpb.setBounds = s.Bounds;
        jpb.setValue  = s.Value;
        jpb.setString = s.String;
        jpb.setIndeterminate = s.Indeterminate;
        varargout = {};
    case 'get'
        [jpb] = varargin{:};
        s.Bounds = jpb.getBounds;
        s.Value  = jpb.getValue;
        s.String = jpb.getString;
        s.Indeterminate = jpb.isIndeterminate;
        varargout = {s};        
end
end %setgetProgbarData()

% ------------------------------------------------------------------------
function varargout = setgetLabelData(opt,varargin)
% Set/get selected progress bar properties
%   progbarData('set',jProgressBar,propstruct)
%   s = progbarData('get',jProgressBar)
switch opt
    case 'set'
        [jtxt,s] = varargin{:};
        jtxt.setBounds = s.Bounds;
        jtxt.setText = s.Text;
        varargout = {};
    case 'get'
        [jtxt] = varargin{:};
        s.Bounds = jtxt.getBounds;
        s.Text = jtxt.getText;
        varargout = {s};        
end
end %setgetLabelData()

% ------------------------------------------------------------------------
function psnUpdateFcn(~,~,h)
% Evaluate all stored position change callback functions
if h.isvalid
    for j = 1:numel(h.psnCbks)
        feval(h.psnCbks{j});
    end
else
    clear h
end
end

% ------------------------------------------------------------------------
function jframe = getJFrame(hFigHandle)

  % Ensure that hFig is a figure handle...
  hFig = ancestor(hFigHandle,'figure');
  if isempty(hFig)
      error(['Cannot retrieve the figure handle for handle ' num2str(hFigHandle)]);
  end

  jframe = [];
  maxTries = 10;
  while maxTries > 0
      try
          % Get the figure's underlying Java frame
          jf = get(handle(hFig),'javaframe');

          % Get the Java frame's root frame handle
          %jframe = jf.getFigurePanelContainer.getComponent(0).getRootPane.getParent;
          try
              jframe = jf.fFigureClient.getWindow;  % equivalent to above...
          catch %#ok<CTCH>
              jframe = jf.fHG1Client.getWindow;  % equivalent to above...
          end
          if ~isempty(jframe)
              break;
          else
              maxTries = maxTries - 1;
              drawnow; pause(0.1);
          end
      catch %#ok<CTCH>
          maxTries = maxTries - 1;
          drawnow; pause(0.1);
      end
  end
  if isempty(jframe)
      error(['Cannot retrieve the java frame for handle ' num2str(hFigHandle)]);
  end
end %getJFrame()

clear all

% addbutton(string,callback)
% addbuttons(string1,callback1,string2,callback2,...)

%% ---------------- Normal methods ---------------- %
close all

% Nargin == 0
h = FigLocker;
h.fig = gcf;
h.lock
pause(1);
h.unlock
h.delete

% Nargin == 1
h = FigLocker(gcf);
h.fig = gcf;
h.lock
pause(1);
h.unlock
h.delete


%% ---------------- Static methods ---------------- %
close all

% Should not fail
FigLocker.Unlock(gcf);  

% Lock/unlock
FigLocker.Lock(gcf);    
pause(1)
FigLocker.Unlock(gcf);  

% Retrieve
FigLocker.Lock(gcf);
h = FigLocker.GetLocker(gcf);
h.lock;
pause(1);
h.unlock
h.delete


%% ---------------- Child objects ---------------- %
close all

% Progress bar
h = FigLocker.Lock(gcf);
h.lock
jpb = h.addprogressbar;
pause(1)
h.unlock

% Text
h = FigLocker.Lock(gcf);
h.lock
h.settext('Busy...');
pause(0.5)
h.unlock

% Progress bar & text
h.lock
jpb = h.addprogressbar;
jpb.setIndeterminate(false)
jpb.setValue(72);
h.settext('72% Complete');
pause(0.5)
h.unlock

% More convenient methods:
h.lock
h.settext('Indeterminate:')
h.setprogress(inf)
pause(0.5)
h.settext('Value:')
h.setprogress(42)
pause(0.5)
h.settext('Indeterminate:')
h.setprogress(NaN)
pause(0.5)

%% ------------------ Controls ---------------------- %
h = FigLocker.Lock(gcf);
dlg = @(varargin)warndlg('That worked!');
unlk = @(varargin)h.unlock;
h.settext('This is a control test')
h.addbuttons('Stop button',dlg,'Go',dlg,'Open',dlg,'Unlock',unlk)
h.setprogress(42)

%% ---------------- Repeated updates ---------------- %
close all
h = FigLocker.Lock(gcf);
h.lock;
h.settext('Checking repeated updates');
jpb = h.addprogressbar;

for j = 1:10
    if mod(j,2) == 1
        jpb.setIndeterminate(false); 
        jpb.setValue(randi(100));
        h.settext('Refresh')
    else
        h.settext(['Update ' num2str(j)])
        jpb.setIndeterminate(true);
    end
    pause(0.2)
end


%% ---------------- Visiblity handling ---------------- %
close all
h = FigLocker.Lock(gcf);
h.lock;
h.settext('Checking visibility handling...')
jpb = h.addprogressbar;
jpb.setIndeterminate(false)
jpb.setValue(42);
pause(0.8)
set(h.fig,'Visible','off')
pause(0.8)
set(h.fig,'Visible','on')
pause(0.8)



%% ---------------- Position handling ---------------- %
close all

h = FigLocker.Lock(gcf);
pause(0.8)
psn = get(h.fig,'Position');
psn = psn/2;
set(h.fig,'Position',psn);



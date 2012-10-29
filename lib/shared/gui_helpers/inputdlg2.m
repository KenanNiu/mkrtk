function answer = inputdlg2(varargin)
%INPUTDLG2 Create an inputdlg instance which utilises keypresses to confirm
% and quit.
%
% Improves on the solution from Mathworks by instead editing the dialog while building:
% http://www.mathworks.com.au/support/solutions/en/data/1-39UWQT/index.html?product=ML&solution=1-39UWQT

if nargin > 1
    tag = varargin{2};
else
    tag = ' ';
end

% Need access to hidden handles:
shh = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

% Outer scope functions/variables
geth = @()findobj(0,'Type','figure','Tag',tag,'Name',tag);
li = [];
cbk = [];

% Timer adds a listener which triggers when figure is made visible
t = timer('TimerFcn',@set_listener,...
    'Period',0.005,...
    'ExecutionMode','fixedSpacing',...
    'TasksToExecute',Inf,...
    'StartDelay',0);
start(t);

% Instigate the dialog:
answer = inputdlg(varargin{:});


    % -----------------------------------
    function set_listener(tobj,~)
        hdlg = geth();
        if ~isempty(hdlg)
            li = addlistener(hdlg,'Visible','PostSet',@set_kpfs);
            stop(tobj)
            delete(tobj)
        end
        
    end %set_listener()

    % -----------------------------------
    function set_kpfs(~,eventdata)
        hdlg = geth();
        okBtn = findobj(hdlg,'Tag','OK');
        editbox = findobj(hdlg,'Tag','Edit');
        if ~isempty(okBtn) && ~isempty(editbox)
            cbk = get(okBtn,'Callback');
            set(editbox,'Callback',@doEnter)
            delete(li)
            set(0,'ShowHiddenHandles',shh)
        end
    end %set_kpfs()

% -----------------------------------
    function doEnter(obj, evd) 
        h = get(obj,'Parent');
        x = get(h,'CurrentCharacter');
        if unicode2native(x) == 13
            feval( cbk, obj, evd); %doCallback(obj,evd);
        end
    end %doEnter()

end %inputdlg2()
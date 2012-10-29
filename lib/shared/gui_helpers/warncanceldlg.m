function ButtonName = warncanceldlg(inputStr,Title)
%WARNCANCELDLG Create warning dialog that allows the user to continue or cancel. 
%
% WARNCANCELDLG uses QUESTDLG to create a question dialog with two response
% options: Ok or Cancel.


% NARGIN check:
if nargin < 1
    inputStr = 'Default string';
end
if nargin < 2
    Title = 'Default Title';
end
validateattributes(inputStr, {'char','cell'}, {'2d'},mfilename);

% Format string
if ~iscell( inputStr )
    warnStr = mat2cell(inputStr,ones(size(inputStr,1),1));
end

Black = [0 0 0]/255;

% Determine figure position
FigPos    = get(0,'DefaultFigurePosition');
FigPos(3) = 267;
FigPos(4) =  70;
FigPos    = getnicedialoglocation(FigPos, get(0,'DefaultFigureUnits'));

% Create figure
hFig=dialog(                                    ...
  'Visible'         ,'off'                      , ...
  'Name'            ,Title                      , ...
  'Pointer'         ,'arrow'                    , ...
  'Position'        ,FigPos                     , ...
  'KeyPressFcn'     ,@doFigureKeyPress          , ...
  'IntegerHandle'   ,'off'                      , ...
  'WindowStyle'     ,'normal'                   , ...
  'HandleVisibility','callback'                 , ...
  'CloseRequestFcn' ,@doDelete                  , ...
  'Tag'             ,Title                        ...
  );

% Set positions:
DefOffset  =10;

IconWidth  =54;
IconHeight =54;
IconXOffset=DefOffset;
IconYOffset=FigPos(4)-DefOffset-IconHeight;  %#ok
IconCMap=[Black;get(hFig,'Color')];  %#ok

% Button dimensions:
BtnWidth  = 56;
BtnHeight = 22;
BtnYOffset=DefOffset;

% Message text:
MsgTxtXOffset=IconXOffset+IconWidth;
MsgTxtYOffset=DefOffset+BtnYOffset+BtnHeight;
% Calculate current msg text width and height. If negative,
% clamp it to 1 since its going to be recalculated/corrected later
% based on the actual msg string
MsgTxtWidth=max(1, FigPos(3)-DefOffset-MsgTxtXOffset-IconWidth);
MsgTxtHeight=max(1, FigPos(4)-DefOffset-MsgTxtYOffset);

% Make buttons:
BtnXOffset=[MsgTxtXOffset; FigPos(3)-DefOffset-BtnWidth];
NumButtons = 2;
BtnString = {'Ok','Cancel'};
BtnTag = {'Ok','Cancel'};
BtnHandle = cell(2,1);
CBString='uiresume(gcbf)';
for j = 1:NumButtons
    BtnHandle{j} = uicontrol(hFig            , ...
    'Style'              ,'pushbutton', ...
    'Position'           ,[ BtnXOffset(1) BtnYOffset BtnWidth BtnHeight ]           , ...
    'KeyPressFcn'        ,@doControlKeyPress , ...
    'Callback'           ,CBString    , ...
    'String'             ,BtnString{j}, ...
    'HorizontalAlignment','center'    , ...
    'Tag'                ,BtnTag{j}     ...
    );
end

% Make message:
MsgHandle=uicontrol(hFig            , ...
  'Style'              ,'text'         , ...
  'Position'           ,[MsgTxtXOffset MsgTxtYOffset 0.95*MsgTxtWidth MsgTxtHeight ]              , ...
  'String'             ,{' '}          , ...
  'Tag'                ,'Question'     , ...
  'HorizontalAlignment','left'         , ...
  'FontWeight'         ,'bold'         , ...
  'BackgroundColor'    ,get(hFig,'Color')  , ...
  'ForegroundColor'    ,Black            ...
  );

[WrapString,NewMsgTxtPos]=textwrap(MsgHandle,warnStr,75);

AxesHandle=axes('Parent',hFig,'Position',[0 0 1 1],'Visible','off');

texthandle=text( ...  
    'Parent'              ,AxesHandle                      , ...
    'Units'               ,'pixels'                        , ...
    'Color'               ,get(BtnHandle{1},'ForegroundColor')   , ...
    'HorizontalAlignment' ,'left'                          , ...
    'FontName'            ,get(BtnHandle{1},'FontName')    , ...
    'FontSize'            ,get(BtnHandle{1},'FontSize')    , ...
    'VerticalAlignment'   ,'bottom'                        , ...
    'String'              ,WrapString                      , ...
    'Interpreter'         ,'none'                          , ...
    'Tag'                 ,'Question'                        ...
    );

textExtent = get(texthandle, 'Extent');

% (g357851)textExtent and extent from uicontrol are not the same. For
% window, extent from uicontrol is larger than textExtent. But on Mac, it
% is reverse. Pick the max value. 
MsgTxtWidth  = max([MsgTxtWidth NewMsgTxtPos(3)+2 textExtent(3)]);
MsgTxtHeight = max([MsgTxtHeight NewMsgTxtPos(4)+2 textExtent(4)]);

MsgTxtXOffset = IconXOffset+IconWidth+DefOffset;
FigPos(3) = max(NumButtons*(BtnWidth+DefOffset)+DefOffset, ...
  MsgTxtXOffset+MsgTxtWidth+DefOffset);

% Center Vertically around icon
if IconHeight>MsgTxtHeight,
  IconYOffset=BtnYOffset+BtnHeight+DefOffset;
  MsgTxtYOffset=IconYOffset+(IconHeight-MsgTxtHeight)/2;
  FigPos(4)=IconYOffset+IconHeight+DefOffset;
  % center around text
else
  MsgTxtYOffset=BtnYOffset+BtnHeight+DefOffset;
  IconYOffset=MsgTxtYOffset+(MsgTxtHeight-IconHeight)/2;
  FigPos(4)=MsgTxtYOffset+MsgTxtHeight+DefOffset;
end

BtnXOffset=[...
    (FigPos(3)-DefOffset)/2-BtnWidth
    (FigPos(3)+DefOffset)/2 ];

set(hFig ,'Position',getnicedialoglocation(FigPos, get(hFig,'Units')));

BtnPos=cellfun(@(bh)get(bh,'Position'), BtnHandle, 'UniformOutput', false);
BtnPos=cat(1,BtnPos{:});
BtnPos(:,1)=BtnXOffset;
BtnPos=num2cell(BtnPos,2);
cellfun(@(bh,pos)set(bh, 'Position', pos), BtnHandle, BtnPos, 'UniformOutput', false);


delete(MsgHandle);

set(texthandle, 'Position',[MsgTxtXOffset MsgTxtYOffset 0]);


IconAxes = axes(                         ...
  'Parent'      ,hFig                  , ...
  'Units'       ,'Pixels'              , ...
  'Position'    ,[IconXOffset IconYOffset IconWidth IconHeight], ...
  'NextPlot'    ,'replace'             , ...
  'Tag'         ,'IconAxes'              ...
  );

set(hFig,'NextPlot','add')

load dialogicons.mat warnIconData warnIconMap;
IconData = warnIconData;
warnIconMap(256,:) = get(hFig,'Color');
IconCMap = warnIconMap;

Img = image('CData',IconData,'Parent',IconAxes);
set(hFig, 'Colormap', IconCMap);
set(IconAxes, ...
  'Visible','off'           , ...
  'YDir'   ,'reverse'       , ...
  'XLim'   ,get(Img,'XData'), ...
  'YLim'   ,get(Img,'YData')  ...
  );

% make sure we are on screen
movegui(hFig)

set(hFig ,'WindowStyle','modal','Visible','on');
drawnow;

uicontrol(BtnHandle{2})

uiwait(hFig);

if ishghandle(hFig)
    ButtonName=get(get(hFig,'CurrentObject'),'String');
    doDelete;
else
    ButtonName = '';
end

    function doFigureKeyPress(~, evd)
        switch evd.Key
            case {'return','space'}
                uiresume(gcbf);
            case 'escape'
                doDelete
        end
    end

    function doControlKeyPress(~, evd)
        switch(evd.Key)
            case {'return'}
                uiresume(gcbf)
            case 'escape'
                doDelete
        end
    end

    function doDelete(varargin)
        delete(hFig)
    end

end %warncanceldlg()


function figure_size = getnicedialoglocation(figure_size, figure_units)
% adjust the specified figure position to fig nicely over GCBF
% or into the upper 3rd of the screen

%  Copyright 1999-2010 The MathWorks, Inc.
%  $Revision: 1.1.6.5 $

parentHandle = gcbf;
convertData.destinationUnits = figure_units;
if ~isempty(parentHandle)
    % If there is a parent figure
    convertData.hFig = parentHandle;
    convertData.size = get(parentHandle,'Position');
    convertData.sourceUnits = get(parentHandle,'Units');  
    c = []; 
else
    % If there is no parent figure, use the root's data
    % and create a invisible figure as parent
    convertData.hFig = figure('visible','off');
    convertData.size = get(0,'ScreenSize');
    convertData.sourceUnits = get(0,'Units');
    c = onCleanup(@() close(convertData.hFig));
end

% Get the size of the dialog parent in the dialog units
container_size = hgconvertunits(convertData.hFig, convertData.size ,...
    convertData.sourceUnits, convertData.destinationUnits, get(convertData.hFig,'Parent'));

delete(c);

figure_size(1) = container_size(1)  + 1/2*(container_size(3) - figure_size(3));
figure_size(2) = container_size(2)  + 2/3*(container_size(4) - figure_size(4));


end


% The main reason for this function is the changing of the icon to a
% warning icon instead of a question icon.  This is done with a timer which
% runs while the dialog is being created.  As soon as the image is present,
% it gets switched to the warning symbol, then program flow continues as
% normal

% 16-Mar-2012, Joshua Martin

% % NARGIN check:
% if nargin < 1
%     String = 'Default string';
% end
% if nargin < 2
%     Title = 'Default Title';
% end
% 
% warnIconData = [];
% warnIconMap = [];
% load('dialogicons.mat','warnIconData','warnIconMap')
% 
% % Make a timer to change the icon:
% tag = 'ChangeIconTimer';
% t = timer('TimerFcn',@changeIcon,...
%     'Tag',tag,...
%     'Period',0.001,...
%     'ExecutionMode','fixedSpacing',...
%     'TasksToExecute',Inf,...
%     'StartDelay',0);
% start(t)
% 
% % Now create the dialog & let it all happen:
% ButtonName = questdlg(String,Title,'Ok','Cancel','Cancel');
% 
% % Just check:
% t = timerfindall('Tag',tag);
% if ~isempty(t)
%     stop(t);
%     delete(t);
% end
% 
% %---------------------------------------------
% % Timer function to change image:
%     function changeIcon(varargin)
%         % This function gets called during object creation while the
%         % QUESTDLG dialog is being built.  We need to wait until all the
%         % components are present before we change the icon data to a
%         % warning icon.
%         
%         % We could to this as nested conditions, but this way is more
%         % efficient because we only have to find each object once.  We do,
%         % however, need the variables declared in the main function so they
%         % persisit:
%         hdlg = findobj(0,'Tag',Title,'Visible','off');
%         if ~isempty(hdlg)                                   % Wait for dialog to be created    
%             hax = findobj(hdlg,'Tag','IconAxes');
%             if ~isempty(hax)                                % Wait for axes creation
%                 himg = findobj(hax,'Type','image');
%                 if ~isempty(himg)                           % Wait for image display
%                     set(himg,'CData',warnIconData);         % Set the image CData
%                     warnIconMap(256,:) = get(hdlg,'Color'); % Adjust the image map
%                     set(hdlg, 'Colormap', warnIconMap);     % Set the colormap
%                     stop(varargin{1})                       % Stop the timer
%                     delete(varargin{1})                     % Delete timer
%                 end
%             end
%         end
%     end %changeIcon
% 
% 
% end %warncanceldlg



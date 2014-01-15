function loa_out = calcLinesOfAction(loa,models)

loa_out = []; % In case of error or early return

% Collect data
handles = guidata(findall(0,'Type','figure','Name','Registration'));

p = get(handles.PhaseSlider,'Value');   % Current phase
n = numel(models);                      % Number of models


istendon = strcmpi( {models.Tag}, loa.Item );
tendon = models( istendon );
others = models( ~istendon );

% Motion:
[q,x] = models2qx(models,quaternion(NaN),NaN);

% Create a default mid-plane
msp = calculate_midplane(tendon);

% Prepare a new figure:
hf = figure('Name','Tendon Line of Action Definition');
ha = axes('Parent',hf);
clrs = get(ha,'ColorOrder');
hold(ha,'on')
%keyboard

hl = FigLocker.Lock(hf);
hl.setprogress(inf)
hl.settext('Calculating data');

% Gather all tendon data
X = {};
for j = numel(tendon.LoRes) : -1 : 1
    X{j} = tendon.LoRes(j).traces;
end

% Calculate midlines:
M = curve_midline(X);
S = midline_surface(M);

hl.unlock;


% Draw data
for j = 1:n
    
    % For the tendon, we draw all the traces so the user can take an
    % average
    if istendon(j)
        models(j).LoRes.plot(ha, 'Color','k')
        
    elseif ~isempty( models(j).HiRes ) %  && models(mj).HiRes.Visible
        cld = models(j).HiRes.transform( q(p,j), x(p,j,:) );
        cld.plot(ha, 'Color',clrs(j,:));%, 'LineStyle','none', cprops{:});
    end
    
end


% Draw plane that the user can manipulate:
set( get(ha,'Children'),'HitTest','off' )   % Turn of hit-test on all other children
hp = drawPlane3d(msp);                      % Draw plane
set(hp,'FaceAlpha',0.4)
draggable(hp,'AllowRotate',true)            % Make interactive

% Final configuration:
view(3)
axis(ha,'equal','tight')
grid on
title(ha,'All tendon segmentations are overlaid on the current phase')
%title(ha,'Do [something] to contine, or close the figure to cancel');

set(hf,'KeyPressFcn',@kpf)

waitfor(hf)
%error('Line of action calculations not yet fully implemented at the user level.')

% ------------------------------------------------------------------------
function p = calculate_midplane(bone)
X = cat(1,bone.LoRes.xyz);
C = princomp(X);
p0 = mean(X,1);
p = [p0 C(:,1)' C(:,2)'];


% ------------------------------------------------------------------------
function kpf(fig,keydata)
%KPF gets run if the user presses a key while the target figure is active.
%   

% Only accept "return" or "enter" as action keys
if ~any( strcmpi( {'return','enter'}, lower(keydata.Key) ) )
    return
end

ax = findobj(get(fig,'Children'), 'Type','axes')

title(ax,'Now select two points on the plane to define the line of Action')

%setpointer(fig, 'fleur');

a = ginput(2)


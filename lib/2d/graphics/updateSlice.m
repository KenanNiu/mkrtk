function updateSlice(handles,sj,pj,varargin)
%SLICEUPDATE Do graphics updates when changing to slice sj
%
% Usage:
%   UPDATESLICE(handles)        
%       Use current slice & phase and refresh traces
%
%   UPDATESLICE(handles,sj)     
%       Update slice, leave phase unchanged, and refresh traces
%
%   UPDATESLICE(handles,sj,pj)  
%       Update both slice & phase and refresh traces
%
%   UPDATESLICE(handles,[],[])  
%       Clear all displayed traces
%
%   UPDATESLICE(handles,sj,pj,'silent') 
%       Do a silent update which does not trigger listeners
%
% In each case except the last, the userdata of the displayed image is updated with the
% serial datenumber returned by NOW(). Any listeners which need to listen
% for a graphics update should listen to the image's UserData property.
%
%
%
% See also ANNOTATIONMANAGER and UPDATETRACES

% Get handle to the image object:
hi = findobj(get(handles.axes1,'Children'),'type','image');


if nargin == 3 && isempty([sj pj])
    % Call:
    %   updateSlice(handles,[],[])
    
    updateTraces(handles,sj,pj);    % Remove all displayed traces
    
else
    % Calls:
    %   updateSlice(handles)
    %   updateSlice(handles,sj)
    %   updateSlice(handles,sj,pj)
    
    % Set/get current slice
    if nargin < 2
        sj = current('slice',handles.axes1); % Use current slice
    else
        current('slice',handles.axes1,sj);   % Set current slice
    end
    
    % Set/get current phase
    if nargin < 3
        pj = current('phase',handles.axes1); % Use current phase:
    else
        current('phase',handles.axes1,pj)    % Set current phase
    end
    
    % If slice/phase change was requested, update:
    %   - Image
    %   - Stack annotation
    %   - Image (slice) annotation
    if nargin > 1
        
        % Get image specs:
        IMG = handles.(activeImageField(handles));
        dims = size(IMG);
        ns = dims(end-1);
        np = dims(end);
        I = currentImage(handles,sj,pj);
        
        % Set image
        set(hi,'CData',I)                       % Update current image data
        
        % Update Stack Annotation:
        annotationManager(handles.StackAnnotation,[],sj,ns,pj,np,[],[]);
        
        % Update Image (slice) annotation:
        infoWBMF(handles.figure1)
    end
    
    % Then always update traces:
    updateTraces(handles,sj,pj);
    
end

silent = nargin > 3 && strcmpi('silent',varargin{1});

% Regardless of number of inputs, we want to trigger any listeners. We do
% this last so that all graphics updates are complete.  The exception is
% the case of a silent update: updateSlice(handles,sj,pj,'silent')
if ~silent
    set( hi, 'UserData', now() );
end


% ------------------------------------------------------------------------
function updateTraces(handles,sj,pj)
%UPDATETRACES Update or plot the traces for slice sj
%
% For the special case, sj == [], this function simply clears all the
% displayed traces, call in either of the following ways:
%   updateTraces(handles,[])
%   updateTraces(handles,[],[])
%
% Usage:
%   updateTraces([])                Clear all traces for this slice/phase
%   updateTraces([],[])             Same as above
%   updateTraces(handles,sj,pj)     Display all traces for this slice/phase
%
%
% See also UPDATESLICE


hf = handles.figure1;

% Delete any old overlaid traces:
tracetag = 'traceobj';
traceobj = getappdata(hf,tracetag);     % This may or may not exist
try_delete(traceobj)
setappdata(hf,tracetag,[])              % Clear the handles

% Special case - clear traces & return:
if isempty(sj)
    return
end

% Or if no traces exist - return
if isempty(handles.traces)
    return
end

% Now show traces, if any exist for this frame:
ids = find(...
    ([handles.traces.Slice] == sj) & ([handles.traces.Phase] == pj) );
if ~isempty(ids)
    traceobj = NaN(numel(ids),1);
    for k = 1:numel(ids)
        % Create Trace object
        t = ids(k);
        traceobj(k) = line(handles.traces(t).x,handles.traces(t).y,...
            'Tag',handles.traces(t).Tag,...
            'Color',handles.traces(t).Color,...
            'LineStyle',handles.traces(t).LineStyle,...
            'Parent',handles.axes1);
        % Add context menu
        set(traceobj(k),'UIContextMenu',handles.htmenu)
    end
    setappdata(hf,tracetag,traceobj);
else
    if isappdata(hf,tracetag)
        rmappdata(hf,tracetag)
    end
end



% ------------------------------------------------------------------------
function try_delete(objs)
%TRY_DELETE Try to delete the specified objects
%
% We could use a try/catch block, but this is formally neater:
for j = 1:numel(objs)
    if ishandle(objs(j))    % If the object still exists
        delete(objs(j))
    end
end
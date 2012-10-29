function hax = calcHelicalAxes(hax,models,opt)
% OPT is used by 

q = [models.q];
x = {models.x};

% If registration hasn't yet been performed, q and x will be empty.  Do
% nothing
if isempty(q) || all(cellfun(@isempty,x))
    return
end

% This is temporary ugliness...:
if (nargin == 3) && strcmpi(opt,'raw') && isfield(models,'qraw')
    q = [models.qraw];
    x = {models.xraw};
end

np = max(cellfun(@numel,{models.q}));    % or simply size(q,1)
for j = 1:numel(x)
    if ~isempty(x{j})
        x{j} = reshape(x{j},np,1,3);
    end
end
x = cell2mat(x);

static = [models.HiRes];

% Process every helical axis
for hj = 1:numel(hax)
    
    obj1.Id = find(strcmp(hax(hj).Item1,{models.Tag}));
    obj2.Id = find(strcmp(hax(hj).Item2,{models.Tag}));
    
    obj1.q = q(:,obj1.Id);
    obj2.q = q(:,obj2.Id);
    
    obj1.x = squeeze(x(:,obj1.Id,:));
    obj2.x = squeeze(x(:,obj2.Id,:));
    
   
    % Reserve space:
    axs = NaN(np,3);
    ang = NaN(np,1);
    pt  = NaN(np,3);
    slide = NaN(np,1);
    
    % Process every phase step:
    for pj = 2:np   % ==> upgrate to parfor later
        
        p1 = pj-1;     % phase of point (1)
        p2 = pj;       % phase of point (2)
        
        A0 = static(obj1.Id);
        B0 = static(obj2.Id);
        
        A1t = A0.transform( obj1.q(p1),obj1.x(p1,:) );
        A2t = A0.transform( obj1.q(p2),obj1.x(p2,:) );
        
        B1t = B0.transform( obj2.q(p1),obj2.x(p1,:) );
        B2t = B0.transform( obj2.q(p2),obj2.x(p2,:) );
        
        [axs(pj,:),ang(pj,:),pt(pj,:),slide(pj,:)] = helicalAxis(A1t, A2t, B1t, B2t);
        
    end
    hax(hj).Axis = axs;
    hax(hj).Angle = ang;
    hax(hj).Point = pt;
    hax(hj).Slide = slide;
end


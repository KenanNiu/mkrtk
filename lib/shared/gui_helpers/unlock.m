function unlock
%UNLOCK Unlock any figures that have been locked with LOCKFIG

set(0,'ShowHiddenHandles','on')
hfigs = get(0,'Children');
set(0,'ShowHiddenHandles','off')

for j = 1:numel(hfigs)
    %lockfig(hfigs(j),'unlock')
    FigLocker.Unlock(hfigs(j));
end

function tf = isdeveloper
% Check if user is developer.
%
% If you're poking around in here, feel free to give yourself developer
% privileges.  But consider yourself warned - you're enabling half-broken
% code.  Enjoy!

tf = true;
return


% Collect the current user name:
if isunix
    [~,name] = unix('whoami');
elseif ispc
    [~,name] = system('echo %UserName%');
end

name = strtrim(name);   % deblank & remove line returns

% Developer user names:
devs = {'josh','local-admin','joshmartin','marto'};

% Compare:
if any(strcmpi(devs,name))
    tf = true;
else
    tf = false;
end

function javaSendKey(hf,key,modifier)
%JAVASENDKEY Send key strokes to figure using java robot.
%
% JAVASENDKEY(HF,KEY) sends the key defined by the character in KEY to the
% figure HF.
%
% JAVASENDKEY(HF,KEY,MODIFIER) sends the key KEY to figure HF with
% modifier(s).  MODIFIER can be a string containing one modifier (such as
% 'shift', or 'command') or it can be a cell containing multiple modifiers,
% eg: {'alt','shift'}
%
% Note that if KEY is a capital character, a 'shift' modifier is
% automaticaly applied:
%
% Example 1:
%
%   javaSendKey(gcf,'s','command')  % save figure
%
% Example 2:
%
%   h = figure; 
%   set(h,'WindowKeyPressFcn',@(hf,ev) title(gca,ev.Character))
%   javaSendKey(h,'k','shift')
%   javaSendKey(h,'K')              % Same as above
%   javaSendKey(h,'w','command')    % Mac close window command; 'w' must be lower case
%

% Joshua Martin, May 30, 2012

dt = 0.150; % pause interval

import java.awt.event.*     % key events
robot = java.awt.Robot;

% Get focus:
figure(hf)
drawnow

% Manage modifiers
if nargin < 3
    modifier = {};
elseif ~iscell(modifier)
    modifier = {modifier};
end

% Manage upper case keys
if upper(key) == key % key is upper case
    modifier = [{'shift'}, modifier(:)'];
end

% Now convert to Java keycodes:
jkey = to_jkey(key);
jmodifier = cellfun(@to_jkey,modifier);

% Depress Modifiers (no need to pause):
for j = 1:numel(jmodifier)
    robot.keyPress( jmodifier(j) );
end
    
% Press & release the primary key:
robot.keyPress(jkey)
pause(dt)
robot.keyRelease(jkey)

% Release modifiers
for j = 1:numel(jmodifier)
    robot.keyRelease( jmodifier(j) );
end

    % --------------------------------
    function jkey = to_jkey(keystring)
        % Mappings:
        switch keystring
            case 'command'              % Mac "command" key is...
                keystring = 'meta';     %  "meta" key in java
        end
        % Convert to java key character code:
        jkey = eval(sprintf('KeyEvent.VK_%s',upper(keystring)));
    end

end
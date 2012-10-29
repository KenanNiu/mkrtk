function configureView3(handles)
% CONFIGUREVIEW3 Configure the view.
%
% This function is called when a new set of models are loaded, or when a
% new session is loaded
%
% See also PHASEDISPLAY

hf = handles.figure1;   % Shorthand

configureSlider3d([],handles.Slider3D)  % Update slider limits/steps
figure1_ResizeFcn(hf,[],handles)        % Adjust slider size/position

% Axes properties:
axis(ha,'equal','tight');   % Make axes tight
grid(ha,'on');              % Grid on
set(ha,'Color','none')      % No background colour
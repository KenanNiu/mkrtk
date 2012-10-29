% Class specification & tests for Bone class.
%
% This class (when written) will supersede the "Models" structure.

% Handle class object?
%   -> implications for loading, registering, smoothing, etc


clear all


% Constructors
b = Bone();

% Properties
b.Tag
b.HiRes
b.LoRes

% Dependent properties
b.q     % gets qsmooth if it's not empty, otherwise gets qraw
b.x     % same for xsmooth

% Hidden properties
b.qraw
b.xraw
b.qsmooth
b.xsmooth

% Methods (logical returns)
b.isregistered

% Methods (value returns)

% Methods (object manipulators)
b.smoothpose(fun)
b.clearpose




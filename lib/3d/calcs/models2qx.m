function [q,x] = models2qx(models,qfill,xfill)
%MODELS2QX Convert the q and x fields of MODELS to the (q,x) form required
% here, where q is NP-by-NO, and x is NP-by-NO-by-3
%
% FILL specifies the value with which undefined elements will be filled
%
% This function is called by Solver and phaseDisplay
%
%
assert( isa(qfill,'quaternion') )

if numel(xfill) == 3
    xfill = xfill(:)';
else
    xfill = [1 1 1]*xfill(1);
end

% Convert to quaternion matrix of size num_phases-by-num_models
q = {models.q};   
x = {models.x};
np = max( cellfun(@numel,q) );
if isempty(np)
    q = [];
    x = [];
    return
end

% Fill any empty space with the fill variables:
q( cellfun(@isempty,q) ) = {qfill(ones(np,1))};
x( cellfun(@isempty,x) ) = {xfill(ones(np,1),:)};

% Now for the x cell array, we must convert all cell entries to
% NP-by-1-by-3 for concatenation into a 3D array in the required form;
% the PERMUTE function does the trick:
permOrder = repmat({[1 3 2]},1,numel(models));   % Cell array of permutation orders
x = cellfun(@permute,x,permOrder,'UniformOutput',false); % Shift into 3rd dimension

% Then finally convert to matrix:
q = [q{:}];       
x = [x{:}];     % Collapsing into matrix results in NP-by-NO-by-3








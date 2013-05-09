function C = cellsmash(A,B)
%CELLSMASH Smash cell arrays A & B into the same cell array
% Both A & B can be cell arrays of the same size, or either A or B can be a
% cell array with the other being a string for insertion into every
% element.  The main usecase is prepending a file path onto a list of
% files, as shown below:
%
% Example:
%   C = cellsmash('/full/path/to/file/',{'file1.mat';'file2.mat';'file3.mat'}) 
%
%   C =
%       '/full/path/to/file/file1.mat'
%       '/full/path/to/file/file2.mat'
%       '/full/path/to/file/file3.mat'  
%       
%
% Other examples:
%
%   cellsmash({'hi';'bye';'fly'}, ' Barry')
%   ans = 
%       'hi Barry'
%       'bye Barry'
%       'fly Barry'
%
%   cellsmash('Barry ',{'thought','said','did'})
%   ans
%       'Barry thought'    'Barry said'    'Barry did'
%
%   cellsmash({'Barry';'Jo';'Bill'},{' can'; ' did'; ' will'})
%   ans =  
%       'Barry can'
%       'Jo did'
%       'Bill will'
%   
%
%   C = cellsmash({'hi ','Barry';'bye ','Jo'})
%   ans = 
%           'hi Barry'
%           'bye Jo'


if nargin == 1  % cellsmash(A)
    for j = 1:size(A,1)
        C{j,1} = [A{j,:}];
    end
elseif iscell(A) && iscell(B)
    assert(isequal(numel(A),numel(B)),'If two cell arrays are provided, they must be the same size')
    C = cellfun(@(a,b)[a,b],A,B,'UniformOutput',false);
elseif iscell(A) && ~iscell(B)
    C = cellfun(@(a) [a,B],A,'UniformOutput',false);
elseif ~iscell(A) && iscell(B)
    C = cellfun(@(b) [A,b],B,'UniformOutput',false);
else
    error('Bad inputs to cellsmash()')
end


end %cellsmash()
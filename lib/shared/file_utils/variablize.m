function str = variablize(str)
% VARIABLIZE Convert string into something safe for a variable name
%
% Input can also be a cell array of strings
%
% First replace all non alpha-numeric permitted characters with spaces:


% If input is a cell, concatenate:
if iscell(str)
    str = str(:)';
    str(2,:) = {'_'};
end
str = [str{:}];

rep_idx = ~( isstrprop(str,'alphanum') | str == '_' );

str(rep_idx) = '_';

x = textscan(str,'%s','delimiter','_','MultipleDelimsAsOne',true);
c = x{1}';
str = [sprintf('%s_',c{1:end-1}),c{end}];

% If the name starts with a digit, prepend some string:
isnum = isstrprop(str,'digit');
if isnum(1)
    str = ['x_' str];
end

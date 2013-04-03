function varargout = find_structs_with_tag(varargin)
% FIND_STRUCTS_WITH_TAG Search all input structs for structs which have
% specified tag(s).  Multiptle tags can be specified, and multiple input
% structures of any dimension can be passed in, in any order.
%
% Note: finds from later structs overwrite finds from earlier structs.
%   Example: If two structures A [1x5] and B [2x3] are passed in in this
%   order, and both structures have tags matching an input TAG, then the
%   find results from B will replace those of A.  This is because elements
%   from arbitrary N-by-M structures A and B cannot be reliably
%   concatenated because they may have different field sets.
%
% Usage:
%   [axis_a] = find_structs_with_tag('Axis_a',struct)  
%
%   [line_f,axis_3] = find_structs_with_tag('Line_f','Axis_3',line_structs,axis_structs,other_structs)
%

tag_tf    = cellfun(@(x)isa(x,'char'),varargin);
struct_tf = cellfun(@(x)isa(x,'struct'),varargin);

tags    = varargin(tag_tf);
structs = varargin(struct_tf);

% Configure output cell:
varargout = cell(1,numel(tags));

% Search for every tag:
for j = 1:numel(tags)
    for k = 1:numel(structs)
        % Collect all structs that have this tag
        % ***Note: if more than one of the input structs
        flags = strcmp({structs{k}.Tag},tags{j});
        if any(flags)
            varargout{j} = structs{k}(flags);
        end
    end
end

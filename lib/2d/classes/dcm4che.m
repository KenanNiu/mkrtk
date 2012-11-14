classdef dcm4che
%DCM4CHE class
%
% Properties (SetAccess = protected, GetAccess = private):
%   DicomObj                Object of type org.dcm4che2.data.BasicDicomObject
%
% Constructors:
%   dcm = dcm4che()         Create empty dcm4che object     
%   dcm = dcm4che(fname)    Create dcm4che object from dicom file
%
% Ordinary Methods:
%   dcm.load(fname)         Load dicom file into initialised object
%
%

% Joshua Martin, 20-Apr-2012

properties (SetAccess = private) %(SetAccess = protected, GetAccess = private)
    FileName = [];
    DicomObj = [];
end

methods
    %--------------------------------------------
    function obj = dcm4che( varargin )    % (constructor)
        if ~load_java_libs()
            error('Could not load dcm4che java libraries. The were not found on the matlab path')
        end
        obj.DicomObj = org.dcm4che2.data.BasicDicomObject;
        if nargin >= 1
            obj = load(obj,varargin{1});
        end
    end %dcm4che()
    
    %--------------------------------------------
    function obj = clear(obj)
        obj.FileName = [];
        obj.DicomObj.clear();
    end %clear()
    
    %--------------------------------------------
    function data = get(obj,property)
        if ischar(property)
            data = get_tag_data(obj.DicomObj,property);
        elseif isstruct(property);
            fields = fieldnames(property);
            for j = 1:numel(fields)
                data.(fields{j}) = get_tag_data(obj.DicomObj,fields{j});
            end
        elseif iscell(property)
            fields = property;
            sargs  = fields;
            [sargs{2,:}] = deal({});
            data = struct(sargs{:});
            for j = 1:numel(fields)
                data(1).(fields{j}) = get_tag_data(obj.DicomObj,fields{j});
            end
        end
    end %get() 
    
    %--------------------------------------------
    function I = getImage(obj)
        rows = get_tag_data(obj.DicomObj,'Rows');
        cols = get_tag_data(obj.DicomObj,'Columns');
        pixeldata = get_tag_data(obj.DicomObj,'PixelData');
        I = reshape(pixeldata, rows, cols)';
    end %getImage()
    
    %--------------------------------------------
    function struct = info(obj)
        % Matlab's DICOMINFO function is about 4-5 times faster than our
        % internal java method, so let's just use that:
        struct = dicominfo(obj.FileName);
        %struct = basic_dicom_object2structure( obj.DicomObj );
    end %info()
    
    %--------------------------------------------
    function tf = isprop(obj,property)
        %ISPROP Test to see if PROPERTY is a property of OBJ
        try
            obj.DicomObj.get(org.dcm4che2.data.Tag.(property));
            tf = true;
        catch ME
            if isequal(ME.identifier,'MATLAB:subscripting:classHasNoPropertyOrMethod')
                tf = false;
            else
                rethrow(ME)
            end
        end
    end %isprop
    
    %--------------------------------------------
    function obj = load(obj,fname)
        assert(isa(fname,'char'),'Filename must be a string')
        % Load file data into object using java
        obj.FileName = fname;
        fis = java.io.FileInputStream(fname);
        din = org.dcm4che2.io.DicomInputStream(...
            java.io.BufferedInputStream(fis));
        din.readDicomObject(obj.DicomObj, -1)
        if obj.DicomObj.size == 0
            % Corrupted file - could not read
            error('\nCorrupted file: %s\n',obj.FileName)
        end
    end %load()
    
    
    
    
end %methods

% Static Methods
methods (Static)
    %--------------------------------------------
    function tf = libsloaded()
        %LIBSLOADED Check if java libraries are loaded
        %   This is the basic check to see if the dcm4che class will be
        %   functional, since it depends on these libraries.
        tf = load_java_libs();
    end %libsloaded()
end %methods(Static)
end %classdef
    

% ========================================================================
%   Helper functions

% ------------------------------------------------------------------------
function dinfo = basic_dicom_object2structure(dcmobj)
% Check recursion:
% dstack = dbstack;
% namestack = {dstack.name};
% rlevel = sum(strcmpi(namestack{1},namestack));
% if rlevel > 1
%     fprintf(2,'> > Recursion Level %d\n',rlevel);
%     if rlevel > 3
%         fprintf(2,'  *** Need to check that this will work\n');
%         keyboard
%     end
% end


% Grab the dictionary:
dict = org.dcm4che2.data.ElementDictionary.getDictionary;
%dict = org.dcm4che2.data.UIDDictionary.getDictionary();

% Instantiate a filemetaInfoIterator & get file metainfo:
itObj = dcmobj.fileMetaInfoIterator;
while itObj.hasNext
    object2field()
end % while

% Instantiate a datasetIterator & get data:
itObj = dcmobj.datasetIterator;
while itObj.hasNext
    object2field()
end

    % ---------------------
    function object2field
        obj = itObj.next;       % Get the next object
        %tobj.isEmpty
        %tobj.tag       %its number
        tagname = dict.nameOf(obj.tag).toCharArray';
        fieldname = tagname2fieldname(tagname);
        
        % tobj.getValueAsString(dcmobj.getSpecificCharacterSet,tobj.length)
        if isempty(fieldname)
            return
            % Private tag
            fieldname = get_private_tag(obj);
        end
        dinfo.(fieldname) = get_tag_data(dcmobj,obj.tag);
    end %object2field()

end %basic_dicom_object2structure

% ------------------------------------------------------------------------
function name = get_private_tag(sdeObj)
% Input is a SimpleDicomElement object, output is a string
tagstrs = textscan( sdeObj.toString.toCharArray, '(%4c,%4c)');
name = sprintf('Private_%s_%s',tagstrs{:});
end %get_private_tag()

% ------------------------------------------------------------------------
function fieldname = tagname2fieldname(tagname)
% Return a field name with forced camel case
% First strip any non-printable characters:
S = tagname;
S(S<32) = [];   % 32 is a space; we need them for the next step
S(S>126) = [];
% Then make initial caps:
S = [' ' S];
letters = isletter(S);
nonLetters = ~letters;
initialCapElements = [0 nonLetters] & [letters 0];
S(initialCapElements) = upper(S(initialCapElements));
% Handle special conversions:
expr = {'(s)'};
for j = 1:numel(expr)
    S = strrep(S,expr{j},'');
end
% Then generic removal of all remanining non-alphanumerics:
S(~isstrprop(S,'alphanum')) = [];
fieldname = S;
end %tagname2fieldname()

% ------------------------------------------------------------------------
function data = get_tag_data(dcmobj,tag)
% http://www.dcm4che.org/docs/dcm4che2-apidocs/org/dcm4che2/data/VR.html
%
% TAG can be either the text string of the attribute, or the numeric ID.

% Check recursion:
% dstack = dbstack;
% namestack = {dstack.name};
% rlevel = sum(strcmpi(namestack{1},namestack));
% if rlevel > 1
%     fprintf(2,'> > Recursion Level %d\n',rlevel);
%     if rlevel > 3
%         fprintf(2,'  *** Need to check that this will work\n');
%         keyboard
%     end
% end


% 
% all = {'AE','AS','AT','code','CS','DA','DS','DT','FD','FL',...
%     'headerLength','IS','LO','LT','OB','OF','OW','padding','PN',...
%     'SH','SL','SQ','SS','ST','TM','UI','UL','UM','UN_SIEMENS','US','UT'};
% numeric = {};
% alfa    = {};
% unhandled = all;

% Get the Constant Field Value (it's unique numerical tag):
% http://www.dcm4che.org/docs/dcm4che2-apidocs/constant-values.html
if ~isnumeric(tag)
    try
    attrId = org.dcm4che2.data.Tag.(tag);
    catch ME
        data = [];  % tag doesn't exist
        return
    end
else
    attrId = tag;
end

if isempty(dcmobj.get(attrId))
    data = NaN;
    return
end

% Value Representation:
%   VR = char( dcmobj.vrOf(attrId) ); 
% The above is true, but private tags have their VR embedded, which may be
% different to the result produced from the code above. 
% The following method works for standard and private tags
VR = char(dcmobj.get(attrId).vr);

if isequal(VR,'SQ')
    % Retrieve a sequence:
    nestedObj = dcmobj.getNestedDicomObject(attrId);
    if isempty(nestedObj)
        data = [];
    else
        data = basic_dicom_object2structure(nestedObj);
    end
    return
else
    % Get the raw byte data:
    rawData = dcmobj.getBytes(attrId);
end

 
if dcmobj.bigEndian
    endian = 'B';
else
    endian = 'L';
end
swap = need_swap( endian );

% Now convert the data with the correct method according to VR:
switch VR
case  {'AE','AS','CS','DA','DT','LO','LT','SH','ST','TM','UI','UT'}
    data = deblankAndStripNulls(char(rawData'));
    %processedAttr = deblankAndStripNulls(char(rawAttr.Data));
    
case {'AT'}
    keyboard
    % For historical reasons don't transpose AT.
    %processedAttr = dcm_typecast(rawAttr.Data, 'uint16', swap);
    
case {'DS', 'IS'}
    data = sscanf(char(rawData), '%f\\');
    %processedAttr = sscanf(char(rawAttr.Data), '%f\\');
    
case {'FL', 'OF'}
    data = dcm_typecast(rawData', 'single', swap);
    %processedAttr = dcm_typecast(rawAttr.Data, 'single', swap)';
     
case 'FD'
    data = dcm_typecast(rawData, 'double', swap);
    %processedAttr = dcm_typecast(rawAttr.Data, 'double', swap)';
    
case 'OB'
    data = rawData;
    %processedAttr = rawAttr.Data';
    
case {'OW', 'US'}
    data = dcm_typecast(rawData,'uint16',swap);
    %processedAttr = dcm_typecast(rawAttr.Data, 'uint16', swap)';
    
    
case 'PN'
    data = parsePerson(deblankAndStripNulls(char(rawData')));
    %processedAttr = parsePerson(deblankAndStripNulls(char(rawAttr.Data)));
    
case 'SL'
    data = dcm_typecast(rawData,'int32',swap); %produces column vector
    %processedAttr = dcm_typecast(rawAttr.Data, 'int32', swap)';
    
case 'SQ'
    %keyboard
    %processedAttr = parseSequence(rawAttr.Data, dictionary);

case 'SS'
    data = dcm_typecast(rawData','int16',swap);
    %processedAttr = dcm_typecast(rawAttr.Data, 'int16', swap)';
    
case 'UL'
    data = dcm_typecast(rawData, 'uint32', swap);
    %processedAttr = dcm_typecast(rawAttr.Data, 'uint32', swap)';
    
case 'UN'
    % It's possible that the attribute contains a private sequence
    % with implicit VR; in which case the Data field contains the
    % parsed sequence.
    if isequal(dcmobj.get(attrId),'SQ')
        keyboard
    else
        data = get_tag_data(dcmobj,attrId);
    end
%     if (isstruct(rawAttr.Data))
%         processedAttr = parseSequence(rawAttr.Data, dictionary);
%     else
%         processedAttr = rawAttr.Data';
%     end

otherwise
    % PS 3.5-1999 Sec. 6.2 indicates that all unknown VRs can be
    % interpretted as UN.  
    processedAttr = rawAttr.Data';

end

        
% Change empty arrays to 0-by-0.
if isempty(data)
  data = reshape(data, [0 0]);
end

end %get_tag_data


% ------------------------------------------------------------------------
function str = deblankAndStripNulls(str)
%DEBLANKANDDENULL  Deblank a string, treating char(0) as a blank.
if (isempty(str))
    return
end
while (~isempty(str) && (str(end) == 0))
    str(end) = '';
end
str = deblank(str);
end

% ------------------------------------------------------------------------
function byteOrder = machineEndian
persistent endian
if ~isempty(endian)
    byteOrder = endian;
    return
end
[~,~,endian] = computer;
byteOrder = endian;
end

% ------------------------------------------------------------------------
function tf = need_swap(attr_endian)
assert(any(strcmp(attr_endian,{'L','B'})))
if isequal( machineEndian, attr_endian)
    tf = false;
else
    tf = true;
end
end %need_swap()
    
% ------------------------------------------------------------------------
function out = dcm_typecast(data, type, varargin)
% modelled off [matlabroot /toolbox/images/iptformats/private/dcm_typecast.m] 
swap = [true varargin{1}];
if swap(end)
    switch class(data)
        case {'uint8','int8'}
            out = swapbytes(typecast(data, type));
        otherwise
            out = typecast(swapbytes(data), type);
    end
else
    out = typecast(data, type);
end
end %dcm_typecast

% ------------------------------------------------------------------------
function tf = load_java_libs()
% Check if dcm4che is loaded.  
% If not, try loading. If it fails, return false

% Check it java files are already loaded into memory
tf = ~isempty( which('org.dcm4che2.io.DicomInputStream') );
if tf
    return
end

% If java files are not loaded, search the path and load them
plist = matlabpath;
[fldr,plist] = strtok(plist,pathsep);    
while ~isempty(fldr)
    % See if files are in the the folder FLDR:
    if ~isempty( dir( [fldr filesep 'dcm4che-core*'] ) )
        break
    end
    [fldr,plist] = strtok(plist,pathsep);
end %while

if ~isempty(fldr)
    % To add all files in this directory (slow):
    %fset = dir([fldr filesep '*.jar']);
    %flist = {fset(~cell2mat({fset.isdir})).name}';
    %for j = 1:numel(flist)
    %    javaaddpath([fldr filesep flist{j}]);
    %end
    
    % To add only the specific files (faster):
    flist = {...                    % At a guess, these seem to be
        'dcm4che-core-*.jar';       % the files that we require.
        'dcm4che-image-*.jar';      % 
        'dcm4che-imageio-*.jar';    % We've used wildcards for version
        'dcm4che-iod-*.jar';        % compatibility, and it'll also get
        'slf4j-*.jar';              % the two slf4f-xxx hits we need.
        'log4j-*.jar'};             % 
    for j = 1:numel(flist)
        finfo = dir([fldr filesep flist{j}]);
        for k = 1:numel(finfo)
            javaaddpath([fldr filesep finfo(k).name])
        end
    end
    % Now that they're all added, flag as loaded:
    tf = true;
end %if

end %dcm4che_isloaded()


function personName = parsePerson(personString)
%PARSEPERSON  Get the various parts of a person name

% A description and examples of PN values is in PS 3.5-2000 Table 6.2-1.

pnParts = {'FamilyName'
           'GivenName'
           'MiddleName'
           'NamePrefix'
           'NameSuffix'};

if (isempty(personString))
    personName = makePerson(pnParts);
    return
end

people = tokenize(personString, '\\');  % Must quote '\' for calls to STRREAD.

personName = struct([]);

for p = 1:length(people)

    % ASCII, ideographic, and phonetic characters are separated by '='.
    components = tokenize(people{p}, '=');
    
    if (isempty(components))
        personName = makePerson(pnParts);
        return   
    end
        
    
    % Only use ASCII parts.
    
    if (~isempty(components{1}))
        
        % Get the separate parts of the person's name from the component.
        componentParts = tokenize(components{1}, '^');

        % The DICOM standard requires that PN values have five or fewer
        % values separated by "^".  Some vendors produce files with more
        % than these person name parts.
        if (numel(componentParts) <= 5)

            % If there are the correct numbers, put them in separate fields.
            for q = 1:length(componentParts)
                
                personName(p).(pnParts{q}) = componentParts{q};
                
            end
            
        else
            
            % If there are more, just return the whole string.
            personName(p).FamilyName = people{p};
            
        end
        
    else
        
        % Use full string as value if no ASCII is present.
        if (~isempty(components))
            personName(p).FamilyName = people{p};
        end
    
    end
    
end
end


function personStruct = makePerson(pnParts)
%MAKEPERSON  Make an empty struct containing the PN fields.

for p = 1:numel(pnParts)
    personStruct.(pnParts{p}) = '';
end
end

function tokens = tokenize( input_string, delimiters )
%TOKENIZE  Divide a string into tokens.
%   TOKENS = TOKENIZE(STRING, DELIMITERS) divides STRING into tokens
%   using the characters in the string DELIMITERS. The result is stored
%   in a single-column cell array of strings.
%
%   Examples: 
%
%   tokenize('The quick fox jumped',' ') returns {'The'; 'quick'; 'fox'; 'jumped'}.
%
%   tokenize('Ann, Barry, Charlie',' ,') returns {'Ann'; 'Barry'; 'Charlie'}.
%
%   tokenize('George E. Forsyth,Michael A. Malcolm,Cleve B. Moler',',') returns
%   {'George E. Forsyth'; 'Michael A. Malcolm'; 'Cleve B. Moler'}

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2008/05/12 21:32:01 $

if (~isempty(input_string))
    tcell = textscan(input_string,'%s',-1,'delimiter',delimiters);
    tokens = tcell{1};
else
    tokens = {};
end

end
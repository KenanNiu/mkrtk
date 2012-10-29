function check_dcm4che

pathstr = '/Users/josh/Desktop/MGR_009/DICOM/';
files = {...
    'IM_0007'
    'IM_0008'
    'IM_0009'
    'IM_0010'
    'IM_0011'
    'IM_0012'
    'IM_0013'};

for j = 1:numel(files)
    
    dcm = dcm4che([pathstr files{j}]);
    
    tic
    jdinfo = dcm.info;
    toc
    tic
    dinfo = dicominfo([pathstr files{j}]);
    toc
    
    %checkRead(dcm)
    
end



function checkRead(obj)

dinfo = dicominfo( obj.FileName );

fields = fieldnames(dinfo);

% Churn through all the fields and check the data is the same:
for j = 1:numel(fields)
    
    field = fields{j};
    if ~obj.isprop(field)    
        fprintf(2,'Field not found in dcm4che: %s\n',field);
        continue        
    end
    
    % Unresolved tags, ignoring for now:
    if any( strcmpi(field,{'FileMetaInformationGroupLength'}) )
        continue
    end
    
    
    mdata = dinfo.(field);  % matlab result
    try
    jdata = obj.get(field); % java result
    catch
        keyboard, 
    end
    
    
    if isequal(mdata,jdata)
        fprintf('PASSED: %s\n',field)
    elseif isempty(jdata)
        fprintf(' *** EMPTY: %s  (empty field in DCM4CHE)\n',field)
    else
        keyboard
    end
    
    
end %for
    
classdef DICOMDir < handle
%DICOMDIR DICOMDir object representing the heirarchical structure of a DICOMDIR file
%   
% USAGE:
%   DD = DICOMDir(FILE) where FILE is the path to a given DICOMDIR file
%
% DICOMDIR methods:
%   gui - will open a graphical user interface to browse the DICOMDIR
%   exportImgsToDCM - will export images in the DICOMDIR to a new folder
%   deleteImgs - will delete images in the DICOMDIR from disk
%   imagesInSeries - return a cell array of full filenames of images in study
%   getOrigFullFilenames - returns the full file names of images in the DICOMDIR
%   retrieveExtraFields - retrieves additional fields into patients, studies, series or images
%
%   parseDicomdir - (Static) returns DICOMDIR patients, studies, series and images separately
%
%
% SAHMSTUDY properties:
%   filepath    - Full path of this DICOMDIR file
%   filename    - DICOMDIR filename (defaults to "DICOMDIR")
%   patients    - Structure array containing all PATIENT information
%   studies     - Structure array containing all STUDY information
%   series      - Structure array containing all SERIES information
%   images      - Structure array containing all IMAGE information
%   studiesMap  - N-by-2 array linking STUDIES to their corresponding PATIENT index
%   seriesMap   - N-by-3 array linking SERIES  to their corresponding PATIENT/STUDY indices
%   imagesMap   - N-by-3 array linking IMAGES  to their corresponding PATIENT/STUDY/SERIES indices

% Written by Sven Holcombe (Oct 6, 2011)
%
% Modified by Josh Martin (Jul, 2012)
%   - Include extra method imagesInSeries() which outputs a cell array of
%       file names when provided with the [patient study series] row vector
%       from the seriesMap.
%   - Add helper function to handle different file separators
    
    properties
        filepath = ''
        filename = 'DICOMDIR'
        patients = struct()
        patientsMap = zeros(0,1)
        studies = struct()
        studiesMap = zeros(0,2)
        series = struct()
        seriesMap = zeros(0,3)
        images = struct()
        imagesMap = zeros(0,4)
    end
    
    methods
        function this = DICOMDir(fname, varargin)
            %DICOMDir    Create DICOMDIR object from the given file.
            [this.filepath,this.filename] = fileparts(fname);
            dcmStruct = DICOMDir.parseDicomdir(fname, varargin{:});
            propsToCopy = {'patients','patientsMap','studies','studiesMap','series','seriesMap','images','imagesMap'};
            for i = 1:length(propsToCopy)
                this.(propsToCopy{i}) = dcmStruct.(propsToCopy{i});
            end
        end
        
        function retrieveExtraFields(this,varargin)
            % RETRIEVEEXTRAFIELDS loads dicominfo tags from the first image per patient/study/series
            %
            % DD.RETRIEVEEXTRAFIELDS(...,FIELDTYPE,FIELDNAMES) takes in parameter/value pairs of
            % FIELDTYPE and FIELDNAMES. FIELDTYPE is one of 'patientFields', 'studyFields',
            % 'seriesFields', or 'imageFields'. FIELDNAMES is a cell array of field names to load
            % from the original dicom files. For example:
            %
            % DD.retrieveExtraFields('seriesFields',{'SeriesDescription'},'patientFields',{'PatientBirthDate'})
            %  will look at the dicominfo of the first image for each patient and extract the
            %  'PatientBirthDat', then look at the first image in each series and retrieve its
            %  'SeriesDescription' field. 'PatientBirthDate' and 'SeriesDescription' will be added
            %  to DD.patients and DD.series, respectively.
            
            IP = inputParser;
            IP.addParamValue('patientFields',{});   IP.addParamValue('studyFields',{});
            IP.addParamValue('seriesFields',{});    IP.addParamValue('imageFields',{});
            IP.KeepUnmatched = true; IP.parse(varargin{:})
            [patInds, stuInds, serInds, imgInds] = deal([]);
            if ~isempty(IP.Results.patientFields),  patInds = find(all(this.imagesMap(:,2:end)==1,2)); end
            if ~isempty(IP.Results.studyFields),    stuInds = find(all(this.imagesMap(:,3:end)==1,2)); end
            if ~isempty(IP.Results.seriesFields),   serInds = find(this.imagesMap(:,4)==1); end
            if ~isempty(IP.Results.imageFields),    imgInds = (1:length(this.images))'; end
            dcmToLoadInds = unique(cat(1,patInds, stuInds, serInds, imgInds));
            filenamesToLoad = this.getOrigFullFilenames(dcmToLoadInds);
            diCell = cell(size(dcmToLoadInds));
            h = waitbar(0,'Loading dicominfo files...');
            for i = 1:length(dcmToLoadInds)
                waitbar(i/length(dcmToLoadInds),h)
                try         diCell{i} = dicominfo(filenamesToLoad{i});
                catch ME,   warning(ME)
                end
            end
            close(h)
            % Apply the requested fields to this.(patients/studies/series/images)
            IPflds = {'patientFields', 'studyFields', 'seriesFields', 'imageFields'}; 
            thisFlds = {'patients','studies','series','images'};
            indsSet = {patInds, stuInds, serInds, imgInds};
            for i = 1:4
                for fn = 1:length(IP.Results.(IPflds{i}))
                    valCell = cellfun(@(di)di.(IP.Results.(IPflds{i}){fn}), diCell(ismember(dcmToLoadInds,indsSet{i})),'UniformOutput',false,'ErrorHandler',@(a,b)[]);
                    [this.(thisFlds{i}).(IP.Results.(IPflds{i}){fn})] = valCell{:};
                end
            end
        end
        
        function gui = gui(this)
            %DICOMDir.gui  Creates a GUI to view the patients/studies/series/images in this DICOMDIR
            
            gui = createInterface;
            [shownIndsStu,shownIndsSer] = deal(1);
            [patCell, stuCell, serCell] = createDisplayCells();
            updateInterface('PATIENT');
            
            function [patCell, stuCell, serCell] = createDisplayCells()
                % Set header rows and make the cell to populate patient/study/series tables
                % PATIENTS
                set(gui.patListTabH, 'Data',cell(0,5),'RearrangeableColumns','on',...
                    'ColumnName',{'P','First','Init.','Last','PatientID','DOB'},...
                    'ColumnWidth',{14,'Auto',20,'Auto','Auto','Auto'});
                patCell = [num2cell(this.patientsMap) ...
                    arrayfun(@(x)x.PatientName.GivenName,this.patients,'UniformOutput',false,'ErrorHandler',@(a,b)'???') ...
                    arrayfun(@(x)x.PatientName.MiddleName,this.patients,'UniformOutput',false,'ErrorHandler',@(a,b)'') ...
                    arrayfun(@(x)x.PatientName.FamilyName,this.patients,'UniformOutput',false,'ErrorHandler',@(a,b)'') ...
                    arrayfun(@(x)x.PatientID,this.patients,'UniformOutput',false,'ErrorHandler',@(a,b)'') ...
                    arrayfun(@(x)datestr(datenum(x.PatientBirthDate,'yyyymmdd'),23),this.patients,'UniformOutput',false,'ErrorHandler',@(a,b)'') ...
                    ];
                % STUDIES
                set(gui.stuListTabH, 'Data',cell(0,5),'RearrangeableColumns','on',...
                    'ColumnName',{'P','St','Accession','Study Desc.','Study Date'},...
                    'ColumnWidth',{14,14,'Auto','Auto','Auto'});
                stuCell = [num2cell(this.studiesMap) ...
                    arrayfun(@(x)x.AccessionNumber,this.studies,'UniformOutput',false,'ErrorHandler',@(a,b)nan) ...
                    arrayfun(@(x)x.StudyDescription,this.studies,'UniformOutput',false,'ErrorHandler',@(a,b)'') ...
                    arrayfun(@(x)datestr(datenum([x.StudyDate x.StudyTime],'yyyymmddHHMMSS')),this.studies,'UniformOutput',false,'ErrorHandler',@(a,b)'')];
                % SERIES
                set(gui.serListTabH, 'Data',cell(0,7),'RearrangeableColumns','on',...
                    'ColumnName',{'P','St','Se','Mod','ImageType','Series Desc.','Series No.','Num Imgs.'},...
                    'ColumnWidth',{15,15,17,25,'Auto','Auto','Auto','Auto'});
                serCell = [num2cell(this.seriesMap) ...
                    arrayfun(@(x)x.Modality,this.series,'UniformOutput',false,'ErrorHandler',@(a,b)'') ...
                    arrayfun(@(x)x.ImageType,this.series,'UniformOutput',false,'ErrorHandler',@(a,b)'???') ...
                    arrayfun(@(x)x.SeriesDescription,this.series,'UniformOutput',false,'ErrorHandler',@(a,b)'?') ...
                    arrayfun(@(x)x.SeriesNumber,this.series,'UniformOutput',false,'ErrorHandler',@(a,b)'?') ...
                    cellfun(@(x)nnz(ismember(this.imagesMap(:,1:3),x,'rows')), num2cell(this.seriesMap,2),'UniformOutput',false)];
            end
            
            function updateInterface(levelStr, src, evnt)
                if nargin==3, set(src,'UserData',evnt.Indices), end
                
                switch levelStr
                    case 'PATIENT'
                        set(gui.patListTabH,'Data',patCell)
                        updateInterface('STUDY')
                        
                    case 'STUDY'
                        set(gui.stuListTabH,'UserData',[]);
                        selRowCols = get(gui.patListTabH,'UserData');
                        if isempty(selRowCols), updateInterface('NO PATIENTS SELECTED'), return; end
                        shownIndsStu = find(ismember(this.studiesMap(:,1), selRowCols(:,1)));
                        set(gui.stuListTabH,'Data',stuCell(shownIndsStu,:),'RowName',shownIndsStu)
                        updateInterface('SERIES')
                        
                    case 'SERIES'
                        set(gui.serListTabH,'UserData',[]);
                        selRowCols = get(gui.stuListTabH,'UserData');
                        if isempty(selRowCols),  updateInterface('NO STUDIES SELECTED'), return; end
                        patStuToMatch = this.studiesMap(shownIndsStu(selRowCols(:,1)),:);
                        shownIndsSer = find(ismember(this.seriesMap(:,1:2), patStuToMatch,'rows'));
                        set(gui.serListTabH,'Data',serCell(shownIndsSer,:),'RowName',shownIndsSer)
                        updateInterface('IMAGES')
                        
                    case 'IMAGES'
                        selRowCols = get(gui.serListTabH,'UserData');
                        if isempty(selRowCols)
                            updateInterface('NO SERIES SELECTED'), return;
                        end
                        patStuSerToMatch = this.seriesMap(shownIndsSer(selRowCols(:,1)),:);
                        imgIndsToShow = find(ismember(this.imagesMap(:,1:3), patStuSerToMatch,'rows'));
                        numImgs = length(imgIndsToShow);
                        set(gui.iconSliderH,'Min',1,'Max',max(1.1,numImgs),'Value',1,'UserData',imgIndsToShow)
                        updateInterface('IMAGE')
                        set([gui.buttons.expSelH gui.buttons.delSelH],'Enable','on')
                        
                    case 'IMAGE'
                        imgIndsToShow = get(gui.iconSliderH,'UserData');
                        slideInd = round(get(gui.iconSliderH,'Value'));
                        i = imgIndsToShow(slideInd);
                        if isfield(this.images(i), 'IconImageSequence')
                            icon = this.images(i).IconImageSequence.Item_1;
                            im = reshape(icon.PixelData, icon.Rows, icon.Columns)';
                        else
                            im = nan(64,64);
                        end
                        set(gui.iconImH, 'CData', im), caxis auto
                        drawnow
                        
                        % SPECIAL CASES WHERE NOTHING WAS SELECTED
                    case 'NO PATIENTS SELECTED'
                        set(gui.stuListTabH,'Data',{},'UserData',[])
                        updateInterface('NO STUDIES SELECTED')
                    case 'NO STUDIES SELECTED'
                        set(gui.serListTabH,'Data',{},'UserData',[])
                        updateInterface('NO SERIES SELECTED')
                    case 'NO SERIES SELECTED'
                        set(gui.iconImH, 'CData', nan(64,64))
                        set([gui.buttons.expSelH gui.buttons.delSelH],'Enable','off')
                end
                
            end
            
            function gui = createInterface()
                % Create user interface for the application and return a structure of handles
                gui = struct();
                gui.Window = figure( ...
                    'Name', 'DICOMDIR', 'NumberTitle', 'off', 'MenuBar', 'none', 'Toolbar', 'none', ...
                    'HandleVisibility', 'on', 'Position', [200 200 800 400], ...
                    'KeyPressFcn',@(src,evnt)onKeyPress(evnt));%, 'CloseRequestFcn',@(src,evnt)onExit());
                % + File menu and Main Axes
                FileMenu = uimenu( gui.Window, 'Label', 'File' );
                uimenu( FileMenu, 'Label', 'Exit', 'Callback', @(src,evnt)delete(gui.Window), 'Accelerator','W');
                
                parentBox = uiextras.HBox('Parent', gui.Window, 'Spacing', 3);
                listsBox = uiextras.VBoxFlex('Parent', parentBox, 'Spacing', 3);
                actionsBox = uiextras.VBox('Parent', parentBox, 'Spacing', 3);
                patBox = uiextras.HBox('Parent', listsBox, 'Spacing', 3);
                stuBox = uiextras.HBox('Parent', listsBox, 'Spacing', 3);
                serBox = uiextras.HBox('Parent', listsBox, 'Spacing', 3);
                gui.patListTabH = uitable('Parent',patBox,'CellSelectionCallback',@(src,evnt)updateInterface('STUDY',src,evnt));
                gui.stuListTabH = uitable('Parent',stuBox,'CellSelectionCallback',@(src,evnt)updateInterface('SERIES',src,evnt));
                gui.serListTabH = uitable('Parent',serBox,'CellSelectionCallback',@(src,evnt)updateInterface('IMAGES',src,evnt));
                
                buttsBox = uiextras.VButtonBox('Parent', actionsBox, 'Spacing', 3);
                gui.buttons.expSelH = uicontrol( 'Parent', buttsBox, 'String', 'Export Selected', 'Callback', @(src,evnt)processSelectedSeries('EXPORT') );
                gui.buttons.delSelH = uicontrol( 'Parent', buttsBox, 'String', 'Delete Selected', 'Callback', @(src,evnt)processSelectedSeries('DELETE') );
                gui.iconAxH = axes('Parent',actionsBox,'Position',[ 0 0 1 1]);
                gui.iconSliderH = uicontrol('Style','Slider','parent',actionsBox,'Max',3,'Min',1,'Value',1,'Callback',@(src,evnt)updateInterface('IMAGE'));
                
                set(parentBox, 'Sizes', [-1 128]);
                set(listsBox, 'Sizes', [-1 -1 -1]);
                set(actionsBox, 'Sizes', [-1 128 10]);
                gui.iconImH = imshow(nan(64,64),'Parent',gui.iconAxH,'InitialMagnification','fit');
            end
            
            function processSelectedSeries(procType)
                selRowCols = get(gui.serListTabH,'UserData');
                if isempty(selRowCols), return, end
                patStuSerToMatch = this.seriesMap(shownIndsSer(selRowCols(:,1)),:);
                selectedImgInds = find(ismember(this.imagesMap(:,1:3), patStuSerToMatch,'rows'));
                switch(procType)
                    case 'DELETE'
                        this.deleteImgs(selectedImgInds)
                    case 'EXPORT'
                        this.exportImgsToDCM(selectedImgInds)
                end
            end
        end
        
        function exportImgsToDCM(this, inds, outputDir)
            %EXPORTIMGSTODCM exports the selected images to a given output directory
            %
            %DD.EXPORTIMGSTODCM(INDS,OUTPUTDIR) exports DD.images(INDS) to the given OUTPUTDIR
            origFilenames = this.getOrigFullFilenames(inds);
            if nargin<3, outputDir = uigetdir; end
            if ~exist(outputDir,'dir'), mkdir(outputDir), end
            newFilenames = strcat(outputDir, filesep, createImageStrings(this, inds), '.dcm');
            h = waitbar(0,sprintf('Copying %d files...',length(newFilenames)));
            for f = 1:length(origFilenames)
                waitbar(f/length(origFilenames),h)
                copyfile(origFilenames{f},newFilenames{f})
            end
            if ishandle(h), close(h), end
        end
        
        function deleteImgs(this, inds, force)
            %DELETEIMGS deletes files from disk listed in the DICOMDIR
            %
            %DD.DELETEIMGS(INDS) deletes the files referenced by DD.images(INDS)
            %
            %DD.DELETEIMGS(INDS,FORCE) when FORCE is TRUE, deletes without first warning the user
            
            origFilenames = this.getOrigFullFilenames(inds);
            if nargin<3 || ~force
                resp = questdlg(sprintf('Are you sure you want delete %d files?',numel(origFilenames)),'Warning','Yes','No','Yes');
                if ~strcmp(resp,'Yes'), return; end
            end
            % When deleting files, deleting multiple at once via DOS is quicker than MATLAB delete
            fChunks = unique([0:20:length(inds) length(inds)]);
            h = waitbar(0,'Deleting files...');
            for f = 1:length(fChunks)-1
                waitbar(f/length(fChunks),h,sprintf('Deleting file chunk %d of %d',f,length(fChunks)))
                subInds = (fChunks(f)+1):fChunks(f+1);
                cmd = sprintf('del "%s" %s & ', origFilenames{subInds});
                system(cmd(1:end-3));
            end
            if ishandle(h), close(h), end
        end
        
        function images = imagesInSeries(this, pss)
            % output cell array of file names of full paths of all images
            % in series specified by the 1-by-3 vector: [PATIENT STUDY SERIES]
            if isempty(pss)
                images = {};
            elseif numel(pss)~=3
                error('Must provide [ PATIENT STUDY SERIES ] vector')
            else
                id = ismember(this.seriesMap,pss,'rows');   % Ensure it's a member of the set
                if all(id==0)                               %
                    images = {};
                    return
                end  
                imginds = ismember(this.imagesMap(:,1:3),pss,'rows');
                images = getOrigFullFilenames(this,imginds);
            end
        end
        
        function origFilenames = getOrigFullFilenames(this, inds)
            if isempty(inds)
                origFilenames = {};
            else
                origFilenames = strcat(this.filepath, filesep, {this.images(inds).ReferencedFileID}');
                origFilenames = standardise_filesep(origFilenames);
            end
        end
        
        function newStrings = createImageStrings(this, inds)
            if isempty(inds)
                newStrings = {};
            else
                serNums = cellfun(@num2str,{this.series(this.imagesMap(inds,3)).SeriesNumber}','UniformOutput',false);
                stuDesc = {this.studies(this.imagesMap(inds,2)).StudyDescription}';
                indNos = arrayfun(@(x)sprintf('%06d',x), 1:length(inds),'UniformOutput',false)';
                newStrings = strcat(stuDesc,'_',serNums,'_',indNos);
            end
        end
        
    end
    
    methods (Static = true)
        function dicomDirStruct = parseDicomdir(fname, varargin)
            %PARSEDICOMDIR    Parse a DICOMDIR file into patients, studies, series, and images.
            %
            % dcm = DICOMDir.parseDicomdir(fname) reads the file FNAME, and generates DCM, a
            % structure with fields:
            %  'patients','patientsMap'
            %  'studies', 'studiesMap'
            %  'series',  'seriesMap'
            %  'images',  'imagesMap'
            % Each <x>Map field has the same rows as the <x> field, and gives a mapping of the
            % heirarchy between images, series, studies, and patients. The first column of <x>Map
            % gives the patient index of that item, the second gives the study index, the third
            % gives the series index. For example, if the 'seriesMap' field returns:
            %     [1     1     1
            %      1     1     2
            %      1     1     3
            %      2     1     1
            %      2     1     2
            %      2     1     3
            %      2     2     1
            %      2     2     2
            %      2     2     3],
            % then the DICOMDIR contained two patients. Patient 1 had only 1 study with 3 series.
            % Patient two had two studies (each with 3 series).
            %
            %Example:
            % dcm = DICOMDir.parseDicomdir('c:/folder/DICOMDIR');
            %
            
            % Written 2011-10-05 by Sven Holcombe
            
            try
                dicomHeader = dicominfo(fname);
            catch me
                if strcmp(me.identifier, 'Images:dicominfo:notDICOM')
                    error(me.identifier, '"%s" is not a DICOMDIR file.', fname);
                end
                error(me.identifier, me.message);
            end
            
            % Gather the full DirectoryRecordSequence list, along with the "Type" of each record
            hdrFields = fieldnames(dicomHeader.DirectoryRecordSequence);
            recTypeStrs = cellfun(@(fld)dicomHeader.DirectoryRecordSequence.(fld).DirectoryRecordType, hdrFields,'UniformOutput',false);
            
            % Extract patients, studies, series, and images
            maskSet = [strcmpi('PATIENT',recTypeStrs), strcmpi('STUDY',recTypeStrs), strcmpi('SERIES',recTypeStrs), strcmpi('IMAGE',recTypeStrs)];
            patientsC = cellfun(@(fld)dicomHeader.DirectoryRecordSequence.(fld), hdrFields(maskSet(:,1)),'UniformOutput',false);
            studiesC  = cellfun(@(fld)dicomHeader.DirectoryRecordSequence.(fld), hdrFields(maskSet(:,2)),'UniformOutput',false);
            seriesC   = cellfun(@(fld)dicomHeader.DirectoryRecordSequence.(fld), hdrFields(maskSet(:,3)),'UniformOutput',false);
            imagesC   = cellfun(@(fld)dicomHeader.DirectoryRecordSequence.(fld), hdrFields(maskSet(:,4)),'UniformOutput',false);
            
            % In case of huge DICOMDIR, clear some memory as we gather data. Turn each cell of
            % structures into a uniform struct array. Where a particular element was missing a
            % field, that field will be created but will be empty
            clear dicomHeader
            pats = unifyCellOfStructs(patientsC);   clear patientsC
            studs = unifyCellOfStructs(studiesC);   clear studiesC
            sers = unifyCellOfStructs(seriesC);     clear seriesC
            imgs = unifyCellOfStructs(imagesC);     clear imagesC
            
            % Build a map from each element to its parent element numbers
            pssiNos = zeros(size(maskSet));
            for i = 1:size(maskSet,1)
                col = find(maskSet(i,:),1);
                if isempty(col),continue; end
                pssiNos(i:end,col) = pssiNos(i,col) + 1;
                pssiNos(i:end,(col+1):end) = 0;
            end
            PatientStudySeriesImageNos = pssiNos(maskSet(:,4),:);
            PatientStudySeriesNos = pssiNos(maskSet(:,3),1:3);
            PatientStudyNos = pssiNos(maskSet(:,2),1:2);
            PatientNos = pssiNos(maskSet(:,1),1);
            
            % Package up as output
            dicomDirStruct = struct('patients',pats,'studies',studs,'series',sers,'images',imgs,...
                'patientsMap',PatientNos, 'studiesMap',PatientStudyNos, 'seriesMap', PatientStudySeriesNos, 'imagesMap',PatientStudySeriesImageNos);
            
            function outStruct = unifyCellOfStructs(cellOfStructs, varargin)
                % Subfunction that takes a cell array of structures, then makes one final structure
                % array (one element per cell), with fields matching those in the original cell of
                % structs.
                fieldNameCell = cellfun(@(s)fieldnames(s), cellOfStructs,'UniformOutput',false);
                flatContents = cellfun(@(s)struct2cell(s), cellOfStructs,'UniformOutput',false);
                flatContents = cat(1,flatContents{:});
                strList = cat(1,fieldNameCell{:});
                [~,listI, grpNo] = unique(strList,'first');
                [~,sortI] = sort(listI);
                outCell = cell(numel(cellOfStructs), length(listI));
                rowNos = cell2mat(cellfun(@(f,num)ones(numel(f),1)*num, fieldNameCell, num2cell(1:numel(cellOfStructs))','UniformOutput',false));
                for rowNo = 1:length(rowNos)
                    outCell{rowNos(rowNo), grpNo(rowNo)} = flatContents{rowNo};
                end
                outStruct = cell2struct(outCell(:,sortI), strList(listI(sortI)), 2);
            end
            
        end
    end
end

%%%%%%%%% Helper functions %%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
function C = standardise_filesep(C)
%STANDARDISE_FILESEP Replace non-patform fileseps with current platform
% fileseps.
% Input can be cell array or string
switch filesep
    case '/'
        oldsep = '\';
    case '\'
        oldsep = '/';        
end
C = strrep(C,oldsep,filesep);
end %standardise_filesep()
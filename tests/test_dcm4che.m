function test_dcm4che
% An example of how to use the dcm4che toolkit from matlab for file
% reading or other dicom functions.

% Step 1 is to download the toolkit from www.dcm4che.org. Select the
% dcm4che2 tookit 'bin' archive, and unzip it.

% Then add the java libraries to your path
% checkjava = which('org.dcm4che2.io.DicomInputStream');
% if isempty(checkjava)
%     libpath = '/Users/josh/Downloads/dcm4che-2.0.26/lib/';
%     javaaddpath([libpath 'dcm4che-core-2.0.26.jar']);
%     javaaddpath([libpath 'dcm4che-image-2.0.26.jar']);
%     javaaddpath([libpath 'dcm4che-imageio-2.0.26.jar']);
%     javaaddpath([libpath 'dcm4che-iod-2.0.26.jar']);
%     javaaddpath([libpath 'slf4j-api-1.4.3.jar']);
%     javaaddpath([libpath 'slf4j-api-1.4.3.jar']);
%     javaaddpath([libpath 'slf4j-log4j12-1.4.3.jar']);
%     javaaddpath([libpath 'log4j-1.2.13.jar']);
% end

% % Now create a reader as follows (other examples on dcm4che wiki)
% %testfile = [matlabroot filesep 'toolbox' filesep 'images' toolbox 'imdemos' filesep 'CT-MONO2-16-ankle.dcm'];
%testfile = '/Applications/MATLAB_R2011b.app/toolbox/images/imdemos/CT-MONO2-16-ankle.dcm';
%fis = java.io.FileInputStream(testfile);
%dcm = org.dcm4che2.data.BasicDicomObject;
% din = org.dcm4che2.io.DicomInputStream(...
%     java.io.BufferedInputStream(fis));
% din.readDicomObject(dcm, -1);
% 
% % You can access the fields like this:
% %width = dcm.getInt(org.dcm4che2.data.Tag.Width)
% rows = dcm.getInt(org.dcm4che2.data.Tag.Rows);
% cols = dcm.getInt(org.dcm4che2.data.Tag.Columns);
% pixeldata = dcm.getInts(org.dcm4che2.data.Tag.PixelData);
% 
% img = reshape(pixeldata, rows,cols);
% figure
% imagesc(img);

Pathname = '/Users/josh/Desktop/MGR_008/';
%Pathname = uigetdir(Pathname,'Select the folder that contains the DICOMDIR file');
[dcmSeries] = loaddcmdir(Pathname);

% Output variables:
Pathname = dcmSeries.Path;
Files = dcmSeries.Images;

drawnow

% Read dicom files:
tic
[X,z,s,dinfo] = readDicoms( Pathname, Files );
toc

return
% Or read with Java
tic
nf = numel(Files);
hw = waitbar(0,'Loading Files');
dcm = org.dcm4che2.data.BasicDicomObject;
for j = 1:nf
    waitbar(j/nf,hw,sprintf('Reading file %d of %d',j,nf));
    fis = java.io.FileInputStream([Pathname Files{j}]);
    din = org.dcm4che2.io.DicomInputStream(...
         java.io.BufferedInputStream(fis));
     din.readDicomObject(dcm, -1);
     rows = dcm.getInt(org.dcm4che2.data.Tag.Rows);
     cols = dcm.getInt(org.dcm4che2.data.Tag.Columns);
     if j == 1
         X = NaN(rows,cols,nf);
     end
     pixeldata = dcm.getInts(org.dcm4che2.data.Tag.PixelData);
     X(:,:,j) = reshape(pixeldata, rows, cols);
end
try close(hw), end
toc

end

%function readDicomsJava( 
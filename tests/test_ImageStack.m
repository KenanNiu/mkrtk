function test_ImageStack

X = ImageStack

testfile = 'MGR_028 ANKLE dynamic tracing.mat';
if exist(testfile)
    d = load(testfile)
    S = d.handles.DICOM;
    
    % Convert from structure into class:
    X = ImageStack(S)
    
end



imgset = cellsmash('/Users/josh/Projects/Tendonopathy/Temp data/',...
    {'Rubber_0001.tif';
    'Rubber_0002.tif';
    'Rubber_0003.tif';
    'Rubber_0004.tif';
    'Rubber_0005.tif';
    'Rubber_0006.tif';
    'Rubber_0007.tif';
    'Rubber_0008.tif'});

profile on
X = ImageStack.LoadFiles(imgset)

profile off
profile viewer

keyboard
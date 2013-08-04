function options = paradigmOpts()

options = [];

options.experimentName = '';
options.modelName = '';
options.stimOnTime = 1;
options.ISI = 3;
options.TR=2.0; % TR used in the scan sequence
options.voxSize = [];
options.voxDims = [];
options.dummyScansEst=[0,9]; % # of dummy scans before and after the trials. Since the prescan dummies are already moved, only show post scan 217-208=9
options.dummyScansVal=[0,7]; % # of dummy scans before and after the trials. Since the prescan dummies are already moved, only show post scan 169-162=7
options.sessCount=5;
options.runCount=10;
options.estRuns=[1,3,5,6,8,10];
options.valRuns=[2,4,7,9];
options.flattenSessions =0;
options.estimationType = 'est';
options.validationType = 'val';
options.estDicomFolderPrefix = 'est_IAPS_';
options.valDicomFolderPrefix = 'val_IAPS_';
options.stimIDMapFilename = 'stimIDMap.mat';
% these two will receive the runType, Session, and run parameters in an sprintf command
options.stimOrderFilenameTemplate = '%sStimOrder_Session%i_Run%i.csv'; 
options.fileListFilenameTemplate = '%sFileList_Session%i_Run%i.csv';
options.paradigmFolder = '';
options.stimFilesFolder = '';
options.featuresFolder = '';
options.regressorsFolder = '';
options.maskFolder = '';
options.outputFolder = '';
options.dicomFolder = '';
options.niftiFolder =  '';

options.doVoxelSelection = false;
options.voxelSelectionMask = '';

% use the simple system to define regressors. This requires one or more
% files that contain headers indicating the regressor names, and assumes
% they are NOT categorical and any interactions are already calculated in
% the spreadsheets
options.useSimpleRegressors = false;

% features options
options.regressorNames = {};
options.regressorInteractions = [];
options.regressorEvents = [];
options.regressorFilenames = {};
options.regressorCategorical = [];
options.regressorCategoricalNames = {};
options.regressorValComp = [];
options.regressorNuissanceEvents = {};

end
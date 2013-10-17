function clearInputFiles(options, subj)

% delete all the paradigm files
paradigmFiles = fullfile(options.paradigmFolder, '*');
delete(paradigmFiles);

% delete the features file
featuresFiles = fullfile(options.featuresFolder, '*');
delete(featuresFiles );

% delete the experiment file
experimentFilename = fullfile(options.outputFolder,options.modelName,subj,sprintf('experiment-%s.mat',date));
delete(experimentFilename);





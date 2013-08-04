function clearInputFiles(options)

% delete all the paradigm files
paradigmFiles = fullfile(options.paradigmFolder, '*');
delete(paradigmFiles);

% delete the features file
featuresFiles = fullfile(options.featuresFolder, '*');
delete(featuresFiles );

% delete the experiment file
experimentFiles = fullfile(options.outputFolder, options.modelName, sprintf('experiment-%s.mat',date));
delete(experimentFiles);





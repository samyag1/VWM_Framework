function [ X ] = readNuissanceFiles(runPath, nuissanceFilenames)

% open the nuissance regressor filenames for this run and make a design
% matrix out of them
X = [];
for curFilenameIdx = 1:numel(nuissanceFilenames)
    
    % get the current nuissance regressor filename
    curFilename = nuissanceFilenames{curFilenameIdx};
    curPath = fullfile(runPath, curFilename);

    % since the filename could be a wildcard for files that have the same
    % prefix, but differ after that per run, see if we can find a single
    % file with the given value
    actualFilenames = dir(curPath);
    actualFilenames = {actualFilenames.name};
    
    for i = 1:numel(actualFilenames)
        
        % create the actual file path
        actualPath = fullfile(runPath, actualFilenames{i});
        
        % load the values from the csv file. Use the spm version because
        % it's smart enough to load movement parameter files, which csvread
        % is too dumb to do.
        curNuissance = spm_load(actualPath);
        
        % concatonate to the design matrix
        X = [X curNuissance];
    end
end

end


function selectBestLambda(dataFolder, xpmt, options)

% load all the correlations (prediction accuracy)
cc = concatSubFields('model','cc',1,dataFolder);

% calculate the mean correlation for each lambda
meanCC = mean(cc,1);

% get the lambda index that was best for the all the voxels
[~,bestLambdaIdx] = max(meanCC);

% determine how many runs were used to estimate these weights
runCount = numel(xpmt.runFold{options.condIdx});

% determine all the mat files in the data folder
dataFilenames = dir(fullfile(dataFolder, '*.mat'));

% now right the weights back out to the estimation files
for curFileIdx = 1:numel(dataFilenames)
    
    % create the current filename
    dataFile = fullfile(dataFolder,dataFilenames(curFileIdx).name);
    
    % read in the model stored in the file, add the weights and write
    % it back out
    load(dataFile, 'model');
    
    % now load all the beta weights, this has the weights estimated for all
    % the values of lambda
    [nFeat, nVox, ~, nCrossValFolds] = size(model.weightsFull);
    nNuissanceFeat = size(model.nuissanceWeightsFull,1);
    
    % create the final weights matrix using the lambdas selected above.
    % Get the matrix of all features by all voxels by all crossvalFolds for
    % the best lambda, then take the mean across cross val folds. Reshape
    % is necessary here instead of squeeze in the case that there's 1
    % feature
    weightsFinal = mean(reshape(model.weightsFull(:,:,bestLambdaIdx,:), nFeat, nVox, nCrossValFolds),3);
    nuissanceWeightsFinal = [];
    if nNuissanceFeat > 0
        nuissanceWeightsFinal = mean(reshape(model.nuissanceWeightsFull(:,:,bestLambdaIdx,:), nNuissanceFeat, nVox, nCrossValFolds),3);
    end
    
    % put the final weights back inito the model
    if options.addBias
        model.weights = single(weightsFinal(1:end-runCount,:));
        model.biasTerms = single(weightsFinal(end-runCount+1:end,:));
    else
        model.weights = single(weightsFinal);
        model.biasTerms = zeros(runCount,nVox);
    end
    if nNuissanceFeat > 0
        model.nuissanceWeights = single(nuissanceWeightsFinal);
    else
        model.nuissanceWeights = [];
    end
    model.lambdasUsed = ones(nVox,1).*bestLambdaIdx;
    
    % for efficiency clear out the full weight matrices
    model.weightsFull = [];
    model.nuissanceWeightsFull = [];
    save(dataFile, 'model', '-append');
end

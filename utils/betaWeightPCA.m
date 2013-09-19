function [componentCoefs, voxelLoadings, explainedVar, coefFeatureNames] = betaWeightPCA(niftisFolder, featureNames, binsFIRUse, binsFIRTotal, featureMask)

% determine how many bins and features we're dealing with
binCount = numel(binsFIRTotal);
featuresPerBin = numel(featureNames);
totalFeatureCount = featuresPerBin*binCount;

if notDefined('featureMask')
    featureMaskUse = [];
else
    % if the length of the featureMask isn't the number of features to use then
    % assume it's a list of indices, and not a logical indexing vector, so
    % create a logical indexing vector
    if numel(featureMask) ~= featuresPerBin
        featureMaskUse = false(featuresPerBin,1);
        featureMaskUse(featureMask) = true;
    else
        featureMaskUse = logical(featureMask);
    end
end

% create the filename of the file containing the beta weights
betasFilename = fullfile(niftisFolder, 'betas.mat');
load(betasFilename, 'betas');

% make sure the values passed in match the beta weights matrix
assert(totalFeatureCount == size(betas,1));

% make the mask that only uses the features specified in the mask, for
% only the FIR bins specified in the binsFIRUse
finalFeatureMask = false(totalFeatureCount,1);
coefFeatureNames = {};
for curBinIdx = 1:numel(binsFIRUse)
    
    % get the current bin number
    curBin = binsFIRUse(curBinIdx);
    
    % determine which cunk of feature weights the current Bin is and set
    % all those features to 1 in the mask so they are plotted
    binChunkNo = find(binsFIRTotal == curBin);
    maskStart = 1 + (binChunkNo-1)*featuresPerBin;
    maskEnd = binChunkNo*featuresPerBin;

    % if there a feature mask, then use it to create the mask for the
    % current bin, and to mask out only the feature names being used
    if ~isempty(featureMaskUse)
        finalFeatureMask(maskStart:maskEnd) = featureMaskUse;
        curBinFeatureNames = strcat(featureNames(featureMaskUse'), sprintf('_Bin-%i', curBin));
    else
        finalFeatureMask(maskStart:maskEnd) = true;
        curBinFeatureNames = strcat(featureNames, sprintf('_Bin-%i', curBin));
    end
    
    coefFeatureNames = [coefFeatureNames curBinFeatureNames];
end

% apply the mask to get the actual betas to do the PCA on
betasUse = betas(finalFeatureMask,:);

% do the pca on the beta weights
[componentCoefs, voxelLoadings, ~, ~, explainedVar] = pca(betasUse');

% create a figure that shows the explainable variance as "scree plot"
figure
plot(explainedVar);
xlabel('Component #');
ylabel('% Variance Explained');

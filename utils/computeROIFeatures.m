function [summaryBetas, summaryBetasPercentages] computeROIFeatures(filenameROIMask, estDir, outputDir, featureNames, binsFIR, betaThreshold, binsToPlot)

% determine how many features and bins there are
featureCount = numel(featureNames);
binCount = numel(binsFIR);
assert(featureCount == (size(betas,1)/binCount));

% get all the estimated beta weights
betas = concatSubFields('model','weights',2,estDir);

% load the indices into the original nifti that were modeled
voxIdxs = concatSubFields('model','voxFit',1,estDir);
if numel(voxIdxs) == 0
    voxIdxs = 1:volumeVoxCount;
end

% read in the roi mask
hdrROIMask = spm_vol(filenameROIMask);
dataROIMask = spm_read_vols(hdrROIMask);

% determine the brain space indices of the ROI
indicesROIMask_brain = find(dataROIMask > 0);
indicesCount = numel(indicesROIMask_brain);

% translate the roi indices into beta map space
indicesROIMask_betas = zeros(indicesCount,1);
for curIdx = 1:indicesCount
    indicesROIMask_betas(curIdx) = find(voxIdxs == indicesROIMask_brain(curIdx));
end

% for each feature determine the number of voxels with a beta value above
% the threshold
summaryBetas = zeros(featureCount,binCount);
for curFeature = 1:featureCount
    for curBinIdx = 1:binCount
        curFeatureBin = curFeature + (curBinIdx-1)*featureCount;
        summaryBetas(curFeature,curBinIdx) = sum(betas(curFeatureBin,indicesROIMask_betas) > betaThreshold);
    end
end

% turn the number of voxels above threhold into a percentage
summaryBetasPercentages = float(summaryBetas) ./ indicesCount;

% plot the bins speecified
binsToPlotCount = numel(binsToPlot);
plotRows = 2;
plotColumns = ceil(binsToPlotCount/plotRows);
for curBinIdx = 1:binsToPlotCount
    
    curBin = binsToPlot(curBinIdx);
    
    figure();
    subplot(plotRows,plotColumns,curBinIdx);
    bar(cummaryBetasPercentages(:,curBinIdx));
    set(gca,'XTickLabel',featureNames);
    set(gca,'XTick',1:featureCount);
    xlim([0,featureCount+1]);
    rotateXLabels(gca, 45);
end
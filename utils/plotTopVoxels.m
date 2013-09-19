function plotTopVoxels(estDir, valDir, outputDir, featureNames, brainDims, binsFIR, voxelsToPlot)

% get the variance accounted for in the validation set
ccVal = concatSubFields('model', 'cc', 2, valDir);

% get all the estimated beta weights
betas = concatSubFields('model','weights',2,estDir);

% load the indices into the original nifti that were modeled
voxIdxs = concatSubFields('model','voxFit',1,estDir);
if numel(voxIdxs) == 0
    voxIdxs = 1:volumeVoxCount;
end

% find the top N voxels in terms of correlation coefficient
[ccValSorted, ccValSortedIdxs] = sort(ccVal(:), 'descend');
NaNCount = sum(isnan(ccValSorted));
ccValClean = ccValSorted(NaNCount+1:end);
ccValCleanIdxs = ccValSortedIdxs(NaNCount+1:end);
topCCVal = ccValClean(1:voxelsToPlot);
topCCValIdxs = ccValCleanIdxs(1:voxelsToPlot);
[topX,topY,topZ] = ind2sub(brainDims, voxIdxs(topCCValIdxs));

featureCount = numel(featureNames);

% save out a mat file with the top voxel values
topBetas = betas(:,topCCValIdxs);
topBetasFilename = fullfile(outputDir, 'topBetas.mat');
save(topBetasFilename, 'topBetas', 'featureNames', 'binsFIR', 'topCCValIdxs', 'topX', 'topY', 'topZ');

% print out bar plots of the top voxels
for curTop = 1:voxelsToPlot
    
    for curBinIdx = 1:numel(binsFIR)
        
        % make an invisible bar plot of the event beats for this voxel
        curBin = binsFIR(curBinIdx);
        curBinStart = 1 + (curBinIdx-1)*featureCount;
        curBinEnd = curBinIdx*featureCount;
        curEventBetas = betas(curBinStart:curBinEnd,topCCValIdxs(curTop));
        f = figure('Visible', 'off');
        set(f, 'WindowStyle', 'docked');
        set (f, 'Units', 'normalized', 'Position', [0,0,1,1]); % maximize it to make it BIG
        bar(curEventBetas);
        set(gca,'XTickLabel',featureNames);
        set(gca,'XTick',1:featureCount);
        xlim([0,featureCount+1]);
        rotateXLabels(gca, 45);
        
        
        % create the filename for the current plot and save it
        fileTitle = sprintf('topVoxel%02i_X%i-Y%i-Z%i_cc%3f_bin%i_BetaPlot', curTop, topX(curTop), topY(curTop), topZ(curTop), topCCVal(curTop), curBin);
        curTopFilename = [fileTitle '.png'];
        curTopFilename = fullfile(outputDir, curTopFilename);
        print('-dpng', curTopFilename);
        curTopFilenameFig = [fileTitle '.fig'];
        curTopFilenameFig = fullfile(outputDir, curTopFilenameFig);
        saveas(f, curTopFilenameFig);
    end
end

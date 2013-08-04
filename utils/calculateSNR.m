function [snr] = calculateSNR(Y, stimOrder, binsFIR)
%function [snr] = calculateSNR(Y, stimOrder, binsFIR)

% SAM TODO - make all maps at all the FIR bin offsets
% SAM TODO - implement the explainable variance algorithm Alex implemented

% the structure to contain all the results
snr = [];

% find the unique stim in the stimOrder, not including 0 which indicates no stim.
uniqueStim = unique(stimOrder);
uniqueStim = uniqueStim(uniqueStim~=0);
uniqueStimCount = numel(uniqueStim);

% the columns of Y are the voxels
voxCount = size(Y,2);

% allocate the matrices for the mean and variance values of all stimuli for
% all voxels
snr.meanAll = zeros(uniqueStimCount, voxCount);
snr.varAll = zeros(uniqueStimCount, voxCount);

% iterate through all the unique stim and calculate the mean and variance
% of the repeats of that unique stim and store them in the matrices
for curStimIdx = 1:uniqueStimCount
    
    % get the current stim
    curStim = uniqueStim(curStimIdx);
    
    % get a matrix of the Y data with each row containing data for 1 repeat
    % of the current stim, and the columns all the voxels
    curStimData = Y(stimOrder == curStim,:);
    
    % take the mean and variance of all the cur stim data and store in the
    % matrices
    snr.meanAll(curStimIdx,:) = mean(curStimData,1);
    snr.varAll(curStimIdx,:) = var(curStimData,1,1);
end

% now to calculate the mean and variance maps, we'll take the mean of both
% across all the stimuli.
snr.meanMap = mean(snr.meanAll,1);
snr.varMap = mean(snr.varAll,1);
snr.snrMap = snr.meanMap ./ snr.varMap;
end
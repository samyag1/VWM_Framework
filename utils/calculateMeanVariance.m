function [results, uniqueStimOrder] = calculateMeanVariance(Y,stimOrder, offsetBinsFIR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% the structure to contain all the results
results = [];

% find the unique stim in the stimOrder, not including 0 which indicates no stim.
uniqueStim = unique(stimOrder);
uniqueStim = uniqueStim(uniqueStim~=0);
uniqueStimOrder = sort(uniqueStim);
uniqueStimCount = numel(uniqueStim);

% the columns of Y are the voxels
voxCount = size(Y,2);
dataCount = size(Y,1);

% allocate the matrices for the mean and variance values of all stimuli for
% all voxels
results.meanAll = zeros(uniqueStimCount, voxCount);
results.varAll = zeros(uniqueStimCount, voxCount);

% iterate through all the unique stim and calculate the mean and variance
% of the repeats of that unique stim and store them in the matrices
for curStimIdx = 1:uniqueStimCount
    
    % get the current stim
    curStim = uniqueStimOrder(curStimIdx);
    
    % get a matrix of the Y data with each row containing data for 1 repeat
    % of the current stim, and the columns all the voxels
    offsetYIdxs = find(stimOrder == curStim) + offsetBinsFIR;
    offsetYIdxs = offsetYIdxs(offsetYIdxs<=dataCount);
    curStimData = Y(offsetYIdxs,:);
    
    % take the mean and variance of all the cur stim data and store in the
    % matrices
    results.meanAll(curStimIdx,:) = mean(curStimData,1);
    results.varAll(curStimIdx,:) = var(curStimData,1,1);
end


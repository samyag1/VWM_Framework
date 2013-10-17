function [nullDist, pVals] = permutationTest(Y,X,weights,addBias,excludeValFeatures,permutationIters,testCC)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[nObs, nVox] = size(Y);

nullDist = zeros(permutationIters, nVox);
for curIter = 1:permutationIters
    
    % shuffle the Y Data
    shuffledY = Y(randperm(nObs),:);
    
    % predict the shuffled data
    curModel = mrifForward(shuffledY,X,weights,addBias,excludeValFeatures);
    
    % store the current prediction accuracy values in the null dist
    nullDist(curIter,:) = curModel.cc;
end

% do a one or two tailed test?
pVals = zeros(nVox);
for curVox = 1:nVox
    % find the number of values in the null distribution that are greater
    % than the test value
    nLarger = find(nullDist(:,curVox) > testCC(curVox));
    
    % divide the number of null dist values that are larger than the test
    % value by the size of the null dist to get the p-value
    pVals(curVox) = numel(nLarger) / permutationIters;
end


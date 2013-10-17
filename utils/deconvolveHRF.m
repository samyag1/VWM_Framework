function [ deconvY, uniqueStim ] = deconvolveHRF(Y, stimOrder, hrf);
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[nDataPoints, nVox] = size(Y);

% determine the unique values in the stimOrder
uniqueStim = unique(stimOrder);
uniqueStim(uniqueStim==0) = [];

% VERY IMPORTANT. We'll want the response amplitudes we're estimating here
% to be in the order of the stimuli in the features matrix, which are the 
% indices in the stimOrder vector, since elsewhere that assumption is made
uniqueStim = sort(uniqueStim);

% create the design matrix that has a column for each stim
nStim = numel(uniqueStim);
X = zeros(nDataPoints,nStim);
for curRow = 1:nDataPoints
    curID = find(uniqueStim == stimOrder(curRow));
    if curID ~= 0
        X(curRow,curID) = 1;
    end
end

% iterate over all the voxels and convolve each one with it's unique HRF,
% then do regress to estimate the amplitudes
deconvY = zeros(nStim,nVox);
for curVox = 1:nVox
    
    % convolve the design matrix with the HRF for the current voxel
    XConv = convolveHRF(X, hrf(:,curVox));
    
    % regress our data onto the convolved design matrix to calculate the
    % beta weights which are the response amplitudes for all the stimuli
%    deconvY(:,curVox) = ols(XConv,Y');
    deconvY(:,curVox) = gls(Y(:,curVox),XConv,4);
end

end


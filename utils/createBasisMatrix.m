function [ basisMatrix ] = createBasisMatrix( X, hrfBasis, amplitudes )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% the number of features in the design matrix be the number of response
% amplitudes passed in
assert(size(X,2) == size(amplitudes,1));

nTimePoints = size(X,1);
nBasis = size(hrfBasis,2);

% multiple the X design matrix by the amplitudes which will give us a
% vector of the amplitude of each time point. This assumes the columns of X
% are in the same order as the rows of amplitudes
ampTimeSeries = X*amplitudes;

% iterate through the basis functions and create the basisMatrix by
% convolving the amplitude time series with each of the basis functions
basisMatrix = zeros(nTimePoints,nBasis);
for curBasis = 1:nBasis
    curConv = conv(ampTimeSeries, hrfBasis(:,curBasis));
    basisMatrix(:,curBasis) = curConv(1:nTimePoints);
end


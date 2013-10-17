function [ convX ] = convolveHRF( X, hrf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[nTimePoints, nCol] = size(X);

% iterate over the columns of the design matrix and convolve each one with
% the hrf passed in
convX = zeros(nTimePoints,nCol);
for curCol = 1:nCol
    curConv = conv(X(:,curCol), hrf);
    convX(:,curCol) = curConv(1:nTimePoints);
end

end


function [ output_args ] = findNaNs( data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


for curCol = 1:size(data,2)

    colNaNs = sum(isnan(data(:,curCol)));
    if colNaNs == 0
        fprintf('\nNo NaNs exist in voxel: %i', curCol);
    end
end


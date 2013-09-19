function [ stimFilenames ] = loadStimFilenames( folderName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% load the stimIDMap
stimIDMapFilename = fullfile(folderName, 'stimIDMap.mat');
load(stimIDMapFilename, 'stimIDMap');

[~,order] = sort(cell2mat(stimIDMap.Values));
stimFilenames = stimIDMap.Keys(order);

end


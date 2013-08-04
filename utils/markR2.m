function [ ev ] = markR2(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dm = bsxfun(@minus,data,nanmean(data,2));

dm = reshape(dm,[],size(dm,3));
data = reshape(data,[],size(data,3));

ev = 1-(nanvar(dm))./nanvar(data);

end


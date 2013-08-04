function segIdx = colonWin(lowerBound,winSize,upperBound);
%----------------------------------------------------------------------
%  function segIdx = colonWin(lowerBound,winSize,upperBound);
%----------------------------------------------------------------------
% Creates chunks of indices of length <winSize> between <lowerBound> and
% <upperBound>, with whater remainder for non-divisible ranges being
% allocated to the last chunk.
%
%INPUT
% <lowerBound>: the first index
% <winSize>:    the size of each chunk of indices
% <upperBound>: the last index
%
%OUTPUT:
% <segIdx>:     a cell array of index vectors, each vector being of
%               length winSize, exluding the last, which may include any
%               remaining indices.
%-----------------------------------------------------------------------
% DES

segments = lowerBound:winSize:upperBound;

segIdx = cell(length(segments)-1,1);
for seg = 1:(length(segments) - 1)
	segIdx{seg} = [segments(seg):segments(seg+1)-1];

end
segIdx{end+1} = [segments(end):upperBound];


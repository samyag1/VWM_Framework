function f = combinerawfiles(files,datatype,filedim,ignorebytes,endian,ix)

% function f = combinerawfiles(files,datatype,filedim,ignorebytes,endian,ix)
%
% <files> is like what is passed to matchwildcards.m.
%   the files that are matched should contain
%   raw data all of the same type and dimensions.
% <datatype> (optional) is the data type for the files.
%   if [] or not supplied, do what loadmatrix.m does
%   (based on the first file).
% <filedim> (optional) is the dimensions for a single file, like
%   what can be passed to loadmatrix.  if [] or not supplied,
%   do what loadmatrix.m does (based on the first file).
% <ignorebytes> (optional) is the number of bytes at beginning
%   to ignore.  if [] or not supplied, default to 0.
% <endian> (optional) is the endian encoding, 'l' | 'b' | 'n'
%   if [] or not supplied, default to the output of 
%   defaultendian.m if it exists and, otherwise, 'n'.
% <ix> (optional) is a vector of indices into the dimensions for
%   a single file.  the meaning is that we want only these elements
%   returned.  if [] or not supplied, default to the normal behavior,
%   which is the return all elements (in the original matrix dimensions).
%   if supplied, return A x B where the indices are along A and where
%   the different files are along B.
%
% return a matrix which is the collection of the data
% in each file along the next available dimension.
%
% the returned format is double.

% deal with input
if ~exist('datatype','var') || isempty(datatype)
  datatype = [];
end
if ~exist('filedim','var') || isempty(filedim)
  filedim = [];
end
if ~exist('ignorebytes','var') || isempty(ignorebytes)
  ignorebytes = 0;
end
if ~exist('endian','var') || isempty(endian)
  if isempty(which('defaultendian'))
    endian = 'n';
  else
    endian = defaultendian;
  end
end
if ~exist('ix','var') || isempty(ix)
  ix = [];
end

% calcs
files = matchwildcards(files);
numfiles = length(files);

% init
if (isempty(datatype) || isempty(filedim)) && numfiles>0
%  [a,b,c] = loadmatrix(files{1},datatype,filedim,ignorebytes,endian);
  [data,header,xyz] = readVolume(files{1});
  if isempty(datatype)
    datatype = spm_type(header.dt(1));
%    datatype = b;
  end
  if isempty(filedim)
    filedim = header.dim;
%    filedim = c;
  end
end

% easy case
if isempty(ix)

  % init
  f = zeros([filedim numfiles]);
  numdims = length([filedim numfiles]);
  
  % process each file
  idx = indexall(numdims,numdims,1);
  for p=1:numfiles
    idx{numdims} = p;  % for speed
%    f(idx{:}) = loadmatrix(files{p},datatype,filedim,ignorebytes,endian);
    f(idx{:}) = readVolume(files{p});
  end

% hard case
else
  
  % figure out range and offset
  ss = cell(1,length(filedim));
  [ss{:}] = ind2sub(filedim,ix);
  mn = min(ss{end});
  mx = max(ss{end});
  offset = prod(filedim(1:end-1)) * (mn-1);
  
  % do it
  f = zeros([length(ix) numfiles]);
  for p=1:numfiles
    temp = loadmatrix(files{p},datatype,filedim,ignorebytes,endian,[mn mx]);
    f(:,p) = temp(ix-offset);
  end

end




  %   if mod(p,200)==1
  %     fprintf(1,'p=%d,',p);
  %   end
  % % report
  % fprintf(1,'\ncombinerawfiles: combined %d files.\n',numfiles);

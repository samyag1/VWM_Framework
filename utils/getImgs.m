function [ imgs ] = getImgs( foldername, sort )
%getImgs GET THOSE IMAGES

origdir = pwd;

cd(foldername)

filenames = dir;
filenames = {filenames.name};
filenames = filenames(3:end);

if sort
    filenames = sort(filenames);
end

for i = 1:length(filenames)
    imgs{i} = imread(filenames{i});
end

cd(origdir)

end


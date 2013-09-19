function [ fig ] = imageScatter( data, images, sizeRatio, plotSize )
%imageScatter Scatter plot of images
%   Generates a scatter plot of images. "data" is a two-column matrix
%   containing the desired X and Y locations of the images on the plot.
%   "images" is a cell structure containing all images. "sizeratio" is the
%   ratio of the maximum image dimension to the maximum plot dimension.
%   "plotSize" is the [height width] dimensions of the plot. 

imageCount = numel(images);

maxImDim = sizeRatio*max(plotSize);
[h, w] = size(images{1});
scale = maxImDim / max([h w]);
imagesScaled = cellfun(@(x)double(imresize(x,scale))/255,images,'UniformOutput',false);

figure
plotMat = zeros([plotSize 3])+1;
data(:,1) = data(:,1)-min(data(:,1));
data(:,1) = data(:,1)/max(data(:,1));
data(:,1) = data(:,1)*(plotSize(2)-size(imagesScaled{1},2)-1) + 1;
data(:,2) = data(:,2)-min(data(:,2));
data(:,2) = data(:,2)/max(data(:,2));
data(:,2) = abs(1-data(:,2))*(plotSize(1)-size(imagesScaled{1},1)-1) + 1;
data = floor(data);
for i=1:imageCount 
    
    %# compute XData/YData vectors of each image
    xrange = data(i,1):(data(i,1)+size(imagesScaled{i},1)-1);
    yrange = data(i,2):(data(i,2)+size(imagesScaled{i},2)-1);
    
    try
        plotMat(yrange,xrange,:) = imagesScaled{i};
    catch
        disp(sprintf('Image Number: %i is not the correct dimensions',i));
    end
    
end
    
imshow(plotMat)

end


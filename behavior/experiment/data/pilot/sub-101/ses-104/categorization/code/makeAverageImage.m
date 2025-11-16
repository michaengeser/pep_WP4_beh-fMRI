% function to make an average image out of a set of images (made by ChatGPT)
% provide folder of images as input
% output is the average image

function [averageImage] = makeAverageImage(imageFolder)

% Get a list of all image files in the specified folder
imageFiles = dir(fullfile(imageFolder, '*.tif'));  % Modify the extension if your images have a different format

% Initialize variables
sumImages = zeros(size(imread(fullfile(imageFolder, imageFiles(1).name))));

% Loop through each image and accumulate
for i = 1:numel(imageFiles)
    % Read the image
    currentImage = double(imread(fullfile(imageFolder, imageFiles(i).name)));

    % Accumulate the images
    sumImages = sumImages + currentImage;
end

% Calculate the average image
averageImage = uint8(sumImages / numel(imageFiles));

% Display the average image
figure;
imshow(averageImage);
title('Average Image');

end
clear;

% This script creates a binary mask of each nuc and saves that image.
load('slicesNucs');
load('nucsNum');
path = '/Users/ebahry/Desktop/dataset/';
%path = ('/home/ella/Desktop/dataset/');
files = dir(path);
files = {files().name};
file = files{4};
img = imread([path file]);
imgSize = size(img);

uniqNucsNum = unique(nucsNum);
uniqNucsNum = uniqNucsNum(uniqNucsNum>0);

% init matrix to save x,y,z starting point of each image
% Those values should be eadded to the x, y, z values in the .obj files: 
xyzStart = zeros(750, 3);
spaceX = 200;
spaceY = 400;

for i=1:size(uniqNucsNum)
    
    num = uniqNucsNum(i);
    [row, column] = find(nucsNum == num); 
    
    % Making the mask size by the number of slices needed:
    uniqInColumn = unique(column);
    minInColumn = min(uniqInColumn);
    maxInColumn = max(uniqInColumn);
    numOfSlices = maxInColumn - minInColumn + 1;
    
    % Making the mask just in the right region:
    nuc = slicesNucs{row(1), column(1)};
    minX = nuc(1,1) - spaceX;
    minY = nuc(1,2) - spaceY;
    
    % Saving the start values of x, y, z:
    xyzStart(i, :) = [minX, minY, minInColumn];
    
    maskImg = false(spaceX*2, spaceY*2, numOfSlices); %(imgSize(1), imgSize(2), numOfSlices);
    
    for j=1:size(row,1)
         nuc = slicesNucs{row(j), column(j)};
         polyg = poly2mask(nuc(:,2)-minY, nuc(:,1)-minX, spaceX*2, spaceY*2); %imgSize(1), imgSize(2));
         maskImg(:,:,column(j)-minInColumn+1) = maskImg(:,:,column(j)-minInColumn+1) + polyg;
    end
    
    for j=1:size(maskImg,3)
        imwrite(maskImg(:,:,j), ['../nucsMasks/nuc' num2str(i) '_' num2str(j) '.jpg']);
    end
    
end

save('xyzStart', 'xyzStart');
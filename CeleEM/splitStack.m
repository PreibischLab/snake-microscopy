% Split tiff image

for i = 1:31

    img = imread('../../microscopyImages/z.tif', i);
    
    imwrite(img, ['../../microscopyImages/z' num2str(i) '.tiff']);
end

% Convert all tif to png
inFolder = '../DWing_registered/';
outFolder = '../DWingPNG/';
imsName = dir([inFolder '*.tif']);

for i = 1:size(imsName)
   img = imread([inFolder imsName(i).name]);
   img = uint8(img);
   
   name = imsName(i).name(1:end-4);
   imwrite(img,[outFolder name '.png']);
    
end
path = '/home/ella/Desktop/dataset';
imgFile = '/full_worm_size10_part';
labelFile = '/Labels_full_worm_size10.tif';

labelFiles = dir([path '/Label_size10_*']);

for i=1:80
    img = imread([path imgFile '1.tif'],i);
    img1 = imread([path labelFile],i);
    imshow(img);
    pause(1);
    close all;
end
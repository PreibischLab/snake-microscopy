path = '../nucsMasks/';
pathWrite = '/Users/ebahry/Desktop/EM/EMOBJzCorrected/';
fileName = 'combinedNuc';
fileType = '.tif';
load('xyzStart.mat');

for i=1:721
    
    imfile = ([path fileName num2str(i) fileType]);
    imInfo = imfinfo(imfile);
    zSize = numel(imInfo);
    
    img = false(400,800,zSize);
    for j=1:zSize
        temp = imread(imfile,j);
        temp = imbinarize(temp);
        img(:,:,j) = temp;
    end
       
    [f,v] = isosurface(img,1/2);
    
    v = round(v + xyzStart(i,:));
    
    fileID = fopen([pathWrite 'nuc' num2str(i) '.obj'],'w');
    for j=1:size(v,1)
        fwrite(fileID,['v ' num2str(v(j,1)) ' ' num2str(v(j,2)) ' ' num2str(v(j,3)*2.25) ' 1.0' char(10)]);
    end
    
    for j=1:size(f,1)
        fwrite(fileID,['f ' num2str(f(j,1)) ' ' num2str(f(j,2)) ' ' num2str(f(j,3)) ' 1.0' char(10)]);
    end
    
    fclose(fileID);
    
end
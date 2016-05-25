function templates = createTemplateImages(sigma, templateSize)
    
    mamaFolder = '/Users/ebahry/Dropbox/DWingImages/';
    templateImg = imread([mamaFolder 'DWingPNG/template_affine.png']);
    templateImg = imgaussfilt(templateImg, sigma);
    
    load('landmarks');
    Xs = landmarks(:,1);
    Ys = landmarks(:,2);
    
    templates = nan(templateSize, templateSize, size(Xs,1));
    for i=1:40
        temStartX = Xs(i)-templateSize/2;
        temStartY = Ys(i)-templateSize/2;
        templates(:,:,i) = templateImg(temStartY:temStartY+templateSize-1, temStartX:temStartX+templateSize-1);
    end

end
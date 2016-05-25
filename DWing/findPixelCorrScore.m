
sigma = 3;
templateSize = 64;
mamaFolder = '/Users/ebahry/Dropbox/DWingImages/';
img = imread([mamaFolder 'DWingPNG/template_affine.png']);
img = imgaussfilt(img, sigma);
mask = imread([mamaFolder 'DWingPNG/template_mask.png']);
% nnz(mask)
[rows, columns] = size(img);

start = 5001;
finish = 18000;

[Xs,Ys] = find(mask);
n = size(Xs,1);

corrScore = nan(rows, columns);
K = finish-start;
startTime = tic;
for i=start:finish
    if mod(i+1-start, 100) == 0
        kk = i - start;
        t = toc(startTime);
        pps = kk/t;
        eta = (K-kk)/pps;
        fprintf('Done %i/%i in %.2fs. ~%.2fs remain.\n', kk, K, t, eta);
    end
    temStartX = Xs(i)-templateSize/2;
    temStartY = Ys(i)-templateSize/2;
    pixelTemplateImage = img(temStartY:temStartY+templateSize-1, temStartX:temStartX+templateSize-1);
    if std(double(pixelTemplateImage(:)))
        pixelCorrImg = templateMatching(img, pixelTemplateImage);

        corrScore(Xs(i),Ys(i)) = sum(pixelCorrImg(:)>0.8);
    end
end

figure; imshow(corrScore); impixelinfo;
save(['corrScore_sig' num2str(sigma) '_TS' num2str(templateSize) '_' num2str(finish)], 'corrScore');
% templates = createTemplateImages(templateImg, landmarks, sigma, templateSize)
% 
% Creates a matrix with the neighboring of each landmarks
%  dimensions = templateSize x templateSize x numberOfLandmarks

function templates = createTemplateImages(templateImg, landmarks, sigma, templateSize)
    
%     templateImg = imgaussfilt(templateImg, sigma);
    G = fspecial('gaussian',[5 5],sigma);
    templateImg = imfilter(templateImg,G,'same');
    
    Xs = landmarks(:,1);
    Ys = landmarks(:,2);
    
    templates = nan(templateSize, templateSize, size(landmarks,1));
    for i=1:size(landmarks,1)
        temStartX = Xs(i)-templateSize/2;
        temStartY = Ys(i)-templateSize/2;
        templates(:,:,i) = templateImg( temStartY:temStartY+templateSize-1,...
                                        temStartX:temStartX+templateSize-1  );
    end

end
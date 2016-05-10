function edgeImg = findEdges(img, sigma, thrSobel)

    img = imgaussfilt(img, sigma);
    
    % Sobel Operator on img with threshold
    [magnitudeImg, directionImg] = imgradient(img, 'sobel');
    magnitudeImg(magnitudeImg < thrSobel) = thrSobel;
    edgeImg = -mat2gray(magnitudeImg);
    
end
%% This function finds correlations through a fourier transform 
% AKA template matching

function c =  templateMatching(img, x, y, templateSize)
% function templateMatching(img, pointsTemplates)

%     img = imread('letters2.jpg');
%     if (ndims(img)>2)
%         img = rgb2gray(img);
%     end
%     [rows, cols] = size(img);
%     template = imread('k.jpeg');
%     if (ndims(template)>2)
%         template = rgb2gray(template);
%     end
    
    % Extract template from template image:
    templateImage = '../../DWingPNG/template_affine.png';
    temStartX = x-templateSize/2;
    temStartY = y-templateSize/2;
    template = templateImage(temStartY:temStartY+templateSize-1, temStartX:temStartX+templateSize-1);
    figure; imshow(template,[]); impixelinfo;

    ix = size(img, 2); 
    iy = size(img, 1);
    tx = size(template, 2); % used for bbox placement
    ty = size(template, 1);

    %// Change - Compute the cross power spectrum
    Gi = fft2(img);
    Gt = fft2(template, iy, ix);
    c = real(ifft2((Gi.*conj(Gt))./abs(Gi.*conj(Gt))));
    figure; imshow(c,[]);

end
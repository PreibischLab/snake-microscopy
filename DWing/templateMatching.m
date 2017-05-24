%% This function finds correlations through a fourier transform 
% AKA template matching

function corrScore =  templateMatching(img, template)
    
    % Initialization
    img = rgb2gray(img);
    img = double(img);
    template = double(template);

    %% 2. correlation calculation
    frameMean = conv2(img,ones(size(template))./numel(template),'same');
    templateMean = mean(template(:));
    corrPartI = conv2(img,fliplr(flipud(template-templateMean)),'same')./numel(template);
    corrPartII = frameMean.*sum(template(:)-templateMean);
    stdFrame = sqrt(conv2(img.^2,ones(size(template))./numel(template),'same')-frameMean.^2);
    stdTemplate = std(template(:));
    corrScore = (corrPartI-corrPartII)./(stdFrame.*stdTemplate);

end
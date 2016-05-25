%% initSnake creates the initial snake according to the parameters
% Then calls a function to start snake

function initSnake(img, sigma, alphas, betas, gammas, templateSize, hoodSize)
    
    load('landmarks');
    
    %% Apply Gaussian 
    img = imgaussfilt(img, sigma);
    figure; 
    
    subplot(1,3,2);
    imshow(img,[]); hold on;
    template = imread('../../DWingImages/DWingPNG/template_affine.png');
    imshow(template); hold on;
    scatter(landmarks(:,1), landmarks(:,2), 'filled', 'b');
    pos = get(gca, 'Position');
    pos(1) = 0;
    pos(3) = 0.33;
    set(gca, 'Position', pos)
    
    subplot(1,3,2);
    imshow(img,[]); hold on;
    pos = get(gca, 'Position');
    pos(1) = 0.33;
    pos(3) = 0.33;
    set(gca, 'Position', pos)
    %% Ploting initial snake
    %title('Initial snake:');
    scatter(landmarks(:,1), landmarks(:,2), 'filled', 'r');
    %labels = cellstr(num2str([1:size(landmarks,1)]'));
    %text(landmarks(:,1), landmarks(:,2), labels,'FontSize',20);
    
    subplot(1,3,3);
    imshow(img,[]); hold on;
    
    pos = get(gca, 'Position');
    pos(1) = 0.66;
    pos(3) = 0.33;
    set(gca, 'Position', pos)
    
    
    %% Init gif image to save snake progress
     snakeGif = 'snakeGif.gif';
%     drawInGif(snakeGif, 1);
    figureChildren = get(gca, 'children');
    delete(figureChildren(1:end-1)); 
    
    templates = createTemplateImages(3,32);
    
    %% Starting snake iterations 
    startIterations(img, alphas, betas, gammas, hoodSize, templates, snakeGif);
        
end
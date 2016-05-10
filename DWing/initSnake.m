%% initSnake creates the initial snake according to the parameters
% Then calls a function to start snake

function initSnake(img, sigma, alphas, betas, gammas, templateSize, hoodSize)
    
    load('landmarks');
    
    %% Ploting initial snake
    title('Initial snake:');
    scatter(landmarks(:,1), landmarks(:,2));
    labels = cellstr(num2str([1:size(landmarks,1)]'));
    text(landmarks(:,1), landmarks(:,2), labels,'FontSize',20);
    pause();
    
    %% Apply Gaussian 
    img = imgaussfilt(img, sigma);
    imshow(img,[]);
    
    %% Init gif image to save snake progress
    snakeGif = 'snakeGif.gif';
    drawInGif(snakeGif, 1);
    
    %% Starting snake iterations 
    startIterations(img, alphas, betas, gammas, templateSize, hoodSize);
        
end
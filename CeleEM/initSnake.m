%% initSnake creates the initial snake according to the parameters
% Then calls a function to start snake

function initSnake(img, method, initShape, sigma, thrSobel, alpha, beta, rho, delta, hoodSize, gamma, slice)
    
    figure;
    imagesc(img);
    colormap('gray')
    axis equal tight
    hold on;
    
%     files = dir(['nuc*' num2str(slice) '*']);
%     files = {files(:).name};
%     for i=1:size(files,2)
%         load(files{i})
%         plot(Xs,Ys,'*g');
%     end
    
    f = gcf;
    f.Position = [230 250 500 600];
    
    % If we chose the method dual snake we need to draw 2 initial snakes
    if strcmp(method, 'DUAL')
        nSnakes2Draw = 2;
    % If we chose other snake methods we only need one snake per object
    else
        nSnakes2Draw = 1;
    end
        
    for i = 1 : nSnakes2Draw
        
        % If we have a dual snake, first init state of 1st snake is saved
        if (i == 2) 
            Xs1 = Xs;
            Ys1 = Ys;
        end
         
        % If the user chose to set the snake manually:
        if strcmp(initShape(i), 'MANUAL')
            
            [Xs, Ys] = setManualSnake(img, 6);
            [Xs, Ys] = resampleSnake(img, Xs, Ys);
        
        % If user chose to have a circle as initial snake:
        else
            radius = str2double(initShape(i));
            
            if (nSnakes2Draw == 1) 
                title({'Please click where you want the center of your snake', 'radius is set to:', num2str(radius), ''})
            elseif (i == 2)
                title({'Please click where you want the center of your balloon snake', 'radius is set to:', num2str(radius), ''})
            else
                title({'Please click where you want the center of your shrinking snake', 'radius is set to:', num2str(radius), ''})
            end
            
            % I NEED TO OPTIMIZE THIS VALUE
            % How does our desired #points determined by the radius 
            nPoints = radius * 2;
            
            % Get user input for center
            
            [centerX, centerY] = ginput(1); 
            
            % Init all snake points - Xs and Ys
            Xs = nan(nPoints, 1);
            Ys = nan(nPoints, 1);

            for j= 1:nPoints
                Xs(j) = centerX + floor(radius*cos(j*2*pi/nPoints) + 0.5);
                Ys(j) = centerY + floor(radius*sin(j*2*pi/nPoints) + 0.5);
            end
        end
    end
    
    %% Ploting initial snake/s
    title({'Thank you! Interpolated snake:', ''})
    if (i == 2)
        plot(Xs1, Ys1, 'b', 'LineWidth', 2)
    end
    plot(Xs,Ys, 'g', 'LineWidth', 2);
    
    %% Find Image Edges - Image Energy 
    % Apply Gaussian filter and sobel operator
    edgeImg = findEdges(img, sigma, thrSobel);
    
    %% Init gif image to save snake progress
    snakeGif = 'snakeGif.gif';
    drawInGif(snakeGif, 1);
    
    %% Starting snake iterations (method according to user input)
    if (strcmp(method,  'GREEDY'))
        % need to edit greedy to fit new format
        % snake = greedy(edgeImg, snake, iterations, hoodSize)
        
    % if user chose KASS, BALLOON, or DUAL snake:
    else
        if (strcmp(method, 'KASS'))
            rho = 0;
        end
%         [Xs, Ys] = snakeIterations(edgeImg, Xs, Ys, alpha(1), beta(1), delta, rho, snakeGif);
        inp = 1;
        while inp ~= 0
            [Xs, Ys] = addNormThenViterbi(edgeImg, Xs, Ys, slice);
            slice = slice + 4
            if (slice/100 < 1)
                zero = '0';
            else
                zero = '';
            end
            img = imread(['/Users/ebahry/Desktop/images_em_head/z' zero num2str(slice) '.tif'],'PixelRegion',{[1250 5000],[3500 8800]});
            edgeImg = findEdges(img, sigma, thrSobel);
            figure; imshow(edgeImg,[]);hold on;
            plot(Xs,Ys);
            
            xs1 = Xs;
            ys1 = Ys;
            
            files = dir(['*' num2str(slice) '*']);
            files = {files(:).name};
            for i=1:size(files,2)
                load(files{i})
                plot(Xs,Ys,'*g');
            end
            
            Xs = xs1;
            Ys = ys1;
%             prompt = 'go?';
%             inp = input(prompt);
%             if inp == 2
                [Xs, Ys] = setManualSnake(img, 6);
                [Xs, Ys] = resampleSnake(img, Xs, Ys);
%             end
        end
%         initBigSnake(edgeImg, Xs, Ys);
%         
%         if (strcmp(method, 'DUAL'))
%             rho = 0;
%             [Xs1, Ys1] = snakeIterations(edgeImg, Xs1, Ys1, alpha(2), beta(2), delta, rho, snakeGif);
%             dual2SingleContour(edgeImg, Xs, Ys, Xs1, Ys1);
%         end
    end
end
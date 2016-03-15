% This function gets 3 arguments:
% img - the original image (with some object to detect)
% n - the number of points the user should select to draw the contour 
% numOfPoints - the number of points wanted in the initial snake

% First, the function display the image and asks the user to pick n points
% that roughly discribe the starting location of the snake.
% After the user entered those points the function uses those starting points to
% create more points (in the amount of numOfPoints) 
% Those new points will be the initial snake shape (created by
% interpolation.)

% The function returns the initial snake - the x and y of numOfPoints points.

function snake = setManualSnake(img, n, numOfPoints)
    
    % Empty vectors Xs and Ys will hold the values of the points entered by
    % the user
    Xs(1:n) = NaN;
    Ys(1:n) = NaN;
    
    % Displaying the img to the user
    figure;
    imshow(img);
    % Using gcp to control figure window size - without it title not showing 
    f = gcf();
    f.Position = [230 250 500 600];
    hold on;
    for i=1:n
        title({'Please click on image to roughly set initial snake contour.', 'Disregard distance between points (points will be interpolated).' , 'No need to completely "close" snake.', 'Amount of points remained to be marked:', num2str(n-i+1), ''});
        % Saving the points entered by the user (one point at each iteration) 
        [Xs(i), Ys(i)] = ginput(1); 
        % Ploting the contour created by the user on the image (for indication)
        plot(Xs,Ys,'--o', 'LineWidth', 2);
    end
    
    % Deletes the ploted line 
    children = get(gca, 'children');
    delete(children(1:n));
    
    title({'Thank you! Initial snake will be interpolated.', ''})
    
    %% Create Snake From User Input Contour
    % At the moment just done linearly %%% TO DO %%%
    
    % Get the distance between each snake point to the next (euclDis)
    % and the sum of the whole snake distance (lengthContour)
    [euclDis, lengthContour] = snakeEuclDistance(Xs, Ys, n);
    
    % Creating the snake by adding points at each segment (each segment is defined by the euclidean distance between two *user input points)
    % pos points to our location on the snake while inserting points.
    pos = 1;
    % Going over all the segments
    for i = 1:n-1
        % Calculating the number of points to put on this segment - 
        % Ratio of this to full numOfPoints is matched to ratio of full length and segment lenght 
        numPoints4Seg = round(euclDis(i)*numOfPoints / lengthContour);
        
        % Finding and puting those points on the segment
        newXs = linspace(Xs(i), Xs(i+1), numPoints4Seg);
        newYs = linspace(Ys(i), Ys(i+1), numPoints4Seg);
        
        % Adding points to snake
        snake(1, pos:pos+numPoints4Seg-1) = newXs;
        snake(2, pos:pos+numPoints4Seg-1) = newYs;
        
        pos = pos + numPoints4Seg;
    end
    % The last iteration - add points on snake closing segment
    i = i+1;
    numPoints4Seg = round(euclDis(i)*numOfPoints / lengthContour);
    newXs = linspace(Xs(i), Xs(1), numPoints4Seg);
    newYs = linspace(Ys(i), Ys(1), numPoints4Seg);
    snake(1, pos:pos+numPoints4Seg-1) = newXs;
    snake(2, pos:pos+numPoints4Seg-1) = newYs;
    
end
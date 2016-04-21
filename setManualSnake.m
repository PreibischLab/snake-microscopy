% This function gets the imige:
% img - the original image (with some object to detect)

% First, the function display the image and asks the user to pick n points
% that roughly discribe the starting location of the snake.
% After the user entered those points the function uses those starting points to
% create more points (in the amount of numOfPoints) 
% Those new points will be the initial snake shape (created by
% interpolation.)

% The function returns the initial snake - the x and y of nPoints points.

function [Xs, Ys] = setManualSnake(img)
    
    % Number of snake points for the user to insert:
    n=10;
    % Init empty vectors Xs and Ys will hold the values of the points entered by
    % the user
    Xs(1:n) = NaN;
    Ys(1:n) = NaN;
    
    % Using gcp to control figure window size - without it title not showing 
    
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
    
    %% Create Snake From User Input Contour
    
    [Xs, Ys] = resampleSnake(img, Xs, Ys);
    
end
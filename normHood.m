%% Function - find the current contour, curvature and image values of the entire local hood
%  To find the minimal value in the hood for the point to move to
%  as part of the greedy snake implementation

function [normedHoodCont, normedHoodCur, normedHoodEdge]  = normHood(img, snake, i, hoodSize, averageDist)
    
    % The neigborhod we check is a square around the current snake point.
    % The size of neighborhood is nXn:
    n = hoodSize*2+1;
    
    % x and y are the location of the current snake point:
    x = round(snake(1,i));
    y = round(snake(2,i));

    % Allocating space for the 3 neighborhood matrices
    % Contour matrix - holds the distance from each neighbor point to the next snake point
    normedHoodCont = zeros(n,n);
    % Curvatur matrix - holds how much curve each point in hood allows
    normedHoodCur = zeros(n,n);
    % Edge matrix - holds how much edge is in each point in hood
    % Since we're already in the sobel image we're just taking image values
    normedHoodEdge = img(x-hoodSize:x+hoodSize, y-hoodSize:y+hoodSize);
    
    % Going over all hood points:
    countX = 1;
    for xx=x-hoodSize:x+hoodSize
        countY = 1;
        for yy=y-hoodSize:y+hoodSize
            
            % Finding the contour and curvature values for the current hood point
            currCont = econt(xx, yy, i, snake, averageDist);
            currCur = ecur(xx, yy, i, snake);
            
            % Setting the contour and curvature values for each hood point:
            normedHoodCont(countX, countY) = currCont;
            normedHoodCur(countX, countY) = currCur;
            
            countY = countY+1;
        end
        countX = countX+1;
    end
    
    % Normalizing the matrices
    normedHoodCont = normedHoodCont/max(normedHoodCont(:));
    normedHoodCur = normedHoodCur/max(normedHoodCur(:));
    if (max(normedHoodEdge(:) ~= 0))
        normedHoodEdge = normedHoodEdge/max(normedHoodEdge(:));
    end
end
% Function returns both the total length of the contour and the distance between each point
% Maybe this function is not needed - can be implimented instead like the resampling KASS snake
 
function [euclDis, lengthContour] = snakeEuclDistance(Xs, Ys, n)
    % Create array with the euclidean distance from each (entered by user) point to the next 
    euclDis(1:n) = NaN;
    for i = 1:n-1
        euclDis(i) = sqrt((Xs(i) - Xs(i+1))^2 + (Ys(i) - Ys(i+1))^2);
    end
    % The last iteration - calc distance from last point to first (close circle)
    i = i + 1;
    euclDis(i) = sqrt((Xs(i) - Xs(1))^2 + (Ys(i) - Ys(1))^2);
    
    % Sum up the perimeter of the snake
    lengthContour = sum(euclDis);
end
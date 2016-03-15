%% Evaluating Contour Energy - "The first order differential" %%
function energy = econt(x, y, s, snake, averageDist)

% Computing the distance of current location to next snake point
if (s<size(snake, 2))
    currDist = sqrt((x - snake(1,s+1))^2 + (y - snake(2,s+1))^2);
else
    currDist = sqrt((x - snake(1,1))^2 + (y - snake(2,1))^2);
end

% Computing for the energy for our point:
energy = abs(averageDist - currDist);

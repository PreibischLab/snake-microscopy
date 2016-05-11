%% Evaluating Contour Energy - "The first order differential" %%
function energy = calcContEnergy(x, y, p, Xs, Ys, adjMat, adjDist)
    energy = 0;
    for adjn= 1:3
        if ~isnan(adjMat(p,adjn))
            % Calculating distance between point and adjacent points
            energy = energy + abs(sqrt((x-Xs(adjMat(p,adjn)))^2 + (y-Ys(adjMat(p,adjn)))^2) - adjDist(p, adjn));
        end
    end
    if (isnan(adjMat(p,adjn)))
        energy = energy / 2;
    else
        energy = energy / 3;
    end

end
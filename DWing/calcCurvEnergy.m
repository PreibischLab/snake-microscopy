%% Evaluating Contour Curvatur - "The second order differential"
function curve = calcCurvEnergy(x, y, p, Xs, Ys, adjMat, adjCurve)
     
    adj1X = Xs(adjMat(p,1));
    adj1Y = Ys(adjMat(p,1));
    adj2X = Xs(adjMat(p,2));
    adj2Y = Ys(adjMat(p,2));

    curve = abs((adj1X - 2*x + adj2X)^2 + (adj1Y - 2*y + adj2Y)^2 - adjCurve(p,1));
    
    if ~isnan(adjMat(p,3))
        adj3X = Xs(adjMat(p,3));
        adj3Y = Ys(adjMat(p,3));
        curve = curve + abs((adj1X - 2*x + adj3X)^2 + (adj1Y - 2*y + adj3Y)^2 - adjCurve(p,2));
        curve = curve + abs((adj2X - 2*x + adj3X)^2 + (adj2Y - 2*y + adj3Y)^2 - adjCurve(p,3));
        curve = curve/3;
    end
end
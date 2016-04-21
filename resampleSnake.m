function [Xs, Ys] = resampleSnake(img, Xs, Ys)

    polyg = double(poly2mask(Xs, Ys, size(img,1), size(img,2)));
    [tempX, tempY] = find(polyg,1);
    
    % There is a rare error I need to understand - what's in tempX and tempY when happens
    if ~(tempX > 0)
        tempX
        tempY
        error('resampling error- need to understand this error');
    end
    
    edgePolyg = bwtraceboundary(polyg, [tempX tempY],'N');
    edgePolyg = edgePolyg(1:end-1,:);
    
    % NEED TO IMPROVE THAT!!!
    % How many points do we want our snake to have (compared to size)
    nPoints = size(edgePolyg,1) / 20;

    sample = floor(size(edgePolyg,1)/nPoints);
%     % In order to keep nPoints size constant we need to resample
%     % from a number of points that devides by nPoints
%     points2delete = size(edgePolyg,1) - sample*nPoints;
%     whichPoints2Delete = randsample(size(edgePolyg,1),points2delete);
%     edgePolyg(whichPoints2Delete,:) = [];

    Xs = edgePolyg(1:sample:end, 2);
    Ys = edgePolyg(1:sample:end, 1);
    
end
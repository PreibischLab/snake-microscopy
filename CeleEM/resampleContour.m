% This function resamples the contour.
% It gets the image size, the contour (Xs,Ys) and the space between points (sample)
% and returns the resampled contour.
function [Xs, Ys] = resampleContour(imgSize, Xs, Ys, sample)

    % Creates a binary image of the polygon from its contour
    polyg = double(poly2mask(Ys, Xs, imgSize(1), imgSize(2)));
    
    % Find the contour of the polygon by traveling on it from a specific point.
    [tempX, tempY] = find(polyg,1);    
    if isempty(tempX)
        tempX
    end
    edgePolyg = bwtraceboundary(polyg, [tempX tempY],'N');
    % Delete the last point - because its equal to first point.
    edgePolyg = edgePolyg(1:end-1,:);

    % return the contour:
    Xs = edgePolyg(1:sample:end, 1);
    Ys = edgePolyg(1:sample:end, 2);
    
end
function [Xs, Ys] = resampleSnake4propagate(img, Xs, Ys, sample)

    polyg = double(poly2mask(Ys, Xs, size(img,1), size(img,2)));
    [tempX, tempY] = find(polyg,1);
    
    edgePolyg = bwtraceboundary(polyg, [tempX tempY],'N');
    edgePolyg = edgePolyg(1:end-1,:);

    Xs = edgePolyg(1:sample:end, 1);
    Ys = edgePolyg(1:sample:end, 2);
    
end
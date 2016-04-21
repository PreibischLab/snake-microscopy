function [pointsEdges] = findPointsEdges(edgeImg, Xmat, Ymat)
    
    n = size(Xmat, 1);
    m = size(Xmat, 2);
    
    pointsEdges = nan(n,m);

    for i = 1:n
        for j = 1:m
            
            pointsEdges(i,j) = edgeImg(Ymat(i,j), Xmat(i,j));
            
        end
    end
end
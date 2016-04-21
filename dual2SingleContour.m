function dual2SingleContour1() %edgeImg, Xs, Ys, Xs1, Ys1)

    
    %% Pre Stuff
    load('exampleDualResult.mat');
    img  = imread('n1.tiff');
    img = findEdges(img, 4, 30);
    
    polyg = double(poly2mask(Xs, Ys, size(img,1), size(img,2)));
    [tempX, tempY] = find(polyg,1);
    edgePolyg = bwtraceboundary(polyg, [tempX tempY],'N');
    edgePolyg = edgePolyg(1:end-1,:);
    nPoints = size(edgePolyg,1) / 30;
    sample = floor(size(edgePolyg,1)/nPoints);
    Xs = edgePolyg(1:sample:end, 2);
    Ys = edgePolyg(1:sample:end, 1);
    
    
    polyg = double(poly2mask(Xs1, Ys1, size(img,1), size(img,2)));
    [tempX, tempY] = find(polyg,1);
    edgePolyg = bwtraceboundary(polyg, [tempX tempY],'N');
    edgePolyg = edgePolyg(1:end-1,:);
    nPoints = size(edgePolyg,1) / 30;
    sample = floor(size(edgePolyg,1)/nPoints);
    Xs1 = edgePolyg(1:sample:end, 2);
    Ys1 = edgePolyg(1:sample:end, 1);
    

    % deciding on n and m
    n = size(Xs,1);
    m = 10;
    
        % Deleting points in big snake so 2 snakes will have same size
    points2delete = size(Xs1,1) - n;
	whichPoints2Delete = randsample(size(Xs1,1),points2delete);
	Xs1(whichPoints2Delete) = [];
    Ys1(whichPoints2Delete) = [];
    
    %Creating a line of points between every closest 2 points in 2 snakes
    figure; imshow(img, []); hold on; impixelinfo;
    
    for i= 1:n
       X(i,:) = linspace(Xs(i),Xs1(i),m); 
       Y(i,:) = linspace(Ys(i),Ys1(i),m); 
       
%        plot(X(i,:), Y(i,:), '--*')
    end
    
    
    %% Start finding path of least cost by Viterbi Algorithm 
    
    % Init:
    % The arbitrary starting point will be at a narrow range
    [minDis, i] = min(sqrt((X(:,1) - X(:,end)).^2 + (Y(:,1) - Y(:,end)).^2));
    % Arraigning X and Y accordingly
    Xmat = round([X(i:end,:); X(1:i-1,:)]);
    Ymat = round([Y(i:end,:); Y(1:i-1,:)]);
    
    % Finding the edge (gradient) at each point:
    pointsEdges = findPointsEdges(img, Xmat, Ymat);
    
    % Init saved path
    savedPathCost = nan(m, m, n);
    savedPathLocation = nan(m, m, n);
    
    % Finding the costs of all options in next step:
    E = findCosts(Xmat, Ymat, pointsEdges(2,:), 2);
    % Maybe can add cost of Eext in first step and add to cost
    savedPathLocation(:,:,1) = 5;
    
    % Finding the min cost when first step = 5 (arbitrary first step)
    % In prevBest matrix 1st D is for 2rd step and 2nd D is 3 step
    prevBest = squeeze(E(5,:,:))';
    
    for i = 3:n-1
        E = findCosts(Xmat, Ymat, pointsEdges(i,:), i);
        % So that the matrix will be: 1st D i-1, 2nd D i, 3rd D i+1
        E = permute(E,[1,3,2]);
        for j = 1:m
            E(:,:,j) = E(:,:,j) + prevBest;
        end
        [minCost, minLoc] = min(E,[],1);
        prevBest = squeeze(minCost);
        savedPathLocation(:,:,i-1) = squeeze(minLoc);
        
    end
    
    bestLoc = nan(n,1);
    
    [tempBestCosts, tempBestLocs] = min(prevBest);
    [tempCost, bestLoc(n-1)] = min(tempBestCosts);
    bestLoc(n) = tempBestLocs(bestLoc(n-1));
    for i = n-2:-1:1
        bestLoc(i) = savedPathLocation(bestLoc(i+1),bestLoc(i+2),i);
    end
    
    for i=1:n
       Xs(i) = Xmat(i,bestLoc(i));
       Ys(i) = Ymat(i,bestLoc(i));
    end
    
    plot(Xs, Ys, '--o')
end
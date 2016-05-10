function dual2SingleContour() %edgeImg, Xs, Ys, Xs1, Ys1)

    close all;

%     % Deleting initial snakes from figure 
%     figureChildren = get(gca, 'children');
%     delete(figureChildren(end-2:end-1)); 

    img = imread('../microscopyImages/z=12.png');
    img = imresize(img,0.4);
    load('14.mat');
    Xs1 = cell2mat(Xss(30));
    Ys1 = cell2mat(Yss(30));
    load('4.mat');
    Xs = cell2mat(Xss(50));
    Ys = cell2mat(Yss(50));
    
    % findEdges gets img, sigma and the threshold for after sobel
    img = findEdges(img, 3.5, 30);
    
    % finding equaly distanced points on the small snake
    polyg = double(poly2mask(Xs, Ys, size(img,1), size(img,2)));
    [tempX, tempY] = find(polyg,1);
    edgePolyg = bwtraceboundary(polyg, [tempX tempY],'N');
    edgePolyg = edgePolyg(1:end-1,:);
    nPoints = size(edgePolyg,1) / 20;
    sample = floor(size(edgePolyg,1)/nPoints);
    Xs = edgePolyg(1:sample:end, 2);
    Ys = edgePolyg(1:sample:end, 1);

    % deciding on n and m
    n = size(Xs,1);
    m = 10;
    
    % Finding the line equation between 2 points in the 2 snakes
    % First finding the norm from the small snake
    [NsX, NsY] = addNormToFindLine(Xs,Ys);
    % Calculating the slope of the desired lines
    slopes = NsY ./ NsX;
    % Finding the constant b of the line equation
    Bs = Ys - slopes .* Xs;
 
    % Finding all the points on big snake
    polyg = double(poly2mask(Xs1, Ys1, size(img,1), size(img,2)));
    [tempX, tempY] = find(polyg,1);
    edgePolyg = bwtraceboundary(polyg, [tempX tempY],'N');
    edgePolyg = fliplr(edgePolyg(1:end-1,:));

    % Init big snake points
    Xs1 = nan(n,1);
    Ys1 = nan(n,1);
    
%     % To prove we found the right slopes and bs (and i did!!! this works!!)
%     for i= 1:size(img,1)
%         Yssss(i) = slopes(3)*i + Bs(3);
%     end
%     figure;imshow(img,[]);hold on;
%     plot((1:734),Yssss,'--o');

    % Find the best corresponding point in big snake
    for i= 1:n
       Xs1(i) = 10000;
       Ys1(i) = 10000;

       for j= 1:size(edgePolyg,1)
           if((slopes(i)*edgePolyg(j,1)+Bs(i)-edgePolyg(j,2) < 3) && (slopes(i)*edgePolyg(j,1)+Bs(i)-edgePolyg(j,2) > -3))
               tempX = edgePolyg(j,1);
               tempY = edgePolyg(j,2);   
               if (sqrt((tempX-Xs(i))^2 + (tempY-Ys(i))^2) < sqrt((Xs1(i)-Xs(i))^2 + (Ys1(i)-Ys(i))^2))
                   % 150 is a hack!!
                   if (((i>1) && (sqrt((Xs1(i-1)-Xs(i))^2 + (Ys1(i-1)-Ys(i))^2)<150)) || (i==1))
                       if i>1
                       end
                        Xs1(i) = tempX;
                        Ys1(i) = tempY;
                   end
               end
           end
       end
    end
    
    % In case we didn't find corresponding point in big snake through the norm line
    % find it by closes point
    % This wasn't tested because no privious run had such edge case
    for i=1:n
        if (Xs1(i)==10000)
            for j=1:size(edgePolyg,1)
                if (sqrt((Xs(i)-edgePolyg(j,1))^2 + (Ys(i)-edgePolyg(j,2))^2) < sqrt((Xs(i)-Xs1(i))^2 + (Ys(i)-Ys1(i))^2))
                    Xs1(i) = edgePolyg(j,1);
                    Ys1(i) = edgePolyg(j,2);
                end
            end
        end
    end
    
   [NsX, NsY] = addExpandingForce (Xs1, Ys1);
    
    Xs1 = round(Xs1 + 40 * NsX);
    Ys1 = round(Ys1 + 40 * NsY);
%     figure;imshow(img,[]);hold on;
%     plot(Xs1,Ys1,'--*');
    
    % finding all points between the two snakes:
%     figure; imshow(img,[]); hold on;
    for i= 1:n
       X(i,:) = linspace(Xs(i),Xs1(i),m); 
       Y(i,:) = linspace(Ys(i),Ys1(i),m); 
       
%         plot(X(i,:), Y(i,:), '--*')
    end
    
    
    %% Start finding path of least cost by Viterbi Algorithm 
    % Init:
    % The arbitrary starting point will be at a narrow range
    [minDis, i] = min(sqrt((X(:,1) - X(:,end)).^2 + (Y(:,1) - Y(:,end)).^2));
    % Arraigning X and Y accordingly
    Xmat = round([X(7:end,:); X(1:7-1,:)]);
    Ymat = round([Y(7:end,:); Y(1:7-1,:)]);
    
    % Finding the edge (gradient) at each point:
    pointsEdges = findPointsEdges(img, Xmat, Ymat);
    
    % Init saved path
    savedPathLocation = nan(m, m, n);
    
    % Finding the costs of all options in next step:
    E = findCosts(Xmat, Ymat, pointsEdges(2,:), 2);
    % Maybe can add cost of Eext in first step and add to cost
    savedPathLocation(:,:,1) = 2;
    
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

    %img = imread('../microscopyImages/z=12.png');
    %img = imresize(img,0.4);
    figure;imshow(img,[]);hold on;plot(Xs, Ys, '--*');
    save('39.mat','Xs','Ys');
end
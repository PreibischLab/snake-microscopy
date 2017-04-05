function [Xs, Ys, XsInner, XsOuter, YsInner, YsOuter] = addNormThenViterbi(orgImg, img, nucContour, lambda)
    
    nucContour = unique(nucContour,'rows', 'stable');
    
    Xs = nucContour(:,1);
    Ys = nucContour(:,2);
    
%      figure; imshow(orgImg,[]); hold on; plot(Ys, Xs, '--*');
%      figure; imshow(img,[]); hold on; plot(Ys, Xs, '--*');
    
    % Init n and m
    n = size(Xs,1);
    m = 6;
    % Init norm length
    normLength = 20;
    
    % Find the norm
    [NsX, NsY] = addNormToFindLine(Xs, Ys);
    
    % Add the (multiplication of the) norm on the outside
    XsOuter = round(Xs + normLength * NsX);
    YsOuter = round(Ys + normLength * NsY);
    % Add the (multiplication of the) norm on the inside
    XsInner = round(Xs - normLength * NsX);
    YsInner = round(Ys - normLength * NsY);
    
    % Plot the inner and outer most search points
%     figure;imshow(orgImg,[]);hold on;
%     plot(XsInner,YsInner,'--*');
%     plot(XsOuter,YsOuter,'--*');
    
    % Add m search points on each line
    for i= 1:n
       X(i,:) = linspace(XsInner(i),XsOuter(i),m); 
       Y(i,:) = linspace(YsInner(i),YsOuter(i),m); 
       
%        plot(Y(i,:), X(i,:), '--*');
    end
    
    %% Start finding path of least cost by Viterbi Algorithm 
    
    % Init:
    % Arraigning X and Y accordingly
    Xmat = round(X);
    Ymat = round(Y);
    
    Xs1 = nan;
    Xs2 = nan;
    Ys1 = nan;
    Ys2 = nan;
    
    % Finding the edge (gradient) at each point:
    pointsEdges = findPointsEdges(img, Xmat, Ymat);
    
    % Since we don't know the best starting point, we repeat the process twice
    % Once from the beginning, once from the middle
    for DPround=1:2

        if (DPround == 2)
            pointsEdges = [pointsEdges(floor(n/2+1):end,:);pointsEdges(1:floor(n/2),:)];
            Xmat = [Xmat(floor(n/2+1):end,:);Xmat(1:floor(n/2),:)];
            Ymat = [Ymat(floor(n/2+1):end,:);Ymat(1:floor(n/2),:)];
        end

        % Init saved path
        savedPathLocation = nan(m, m, n);

        % Finding the costs of all options in next step:
        E = findCosts(Xmat, Ymat, pointsEdges(2,:), 2, lambda);
        % Maybe can add cost of Eext in first step and add to cost
        savedPathLocation(:,:,1) = 5;
            

        % Finding the min cost when first step = 5 (arbitrary first step)
        % In prevBest matrix 1st D is for 2rd step and 2nd D is 3 step
        prevBest = squeeze(E(5,:,:))';

        for i = 3:n-1
            E = findCosts(Xmat, Ymat, pointsEdges(i,:), i, lambda);
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
        
        if (DPround==1)
            for i=1:n
               Xs1(i,1) = Xmat(i,bestLoc(i));
               Ys1(i,1) = Ymat(i,bestLoc(i));
            end
        else
            for i=1:n
               Xs2(i,1) = Xmat(i,bestLoc(i));
               Ys2(i,1) = Ymat(i,bestLoc(i));
            end
        end
    end
    
    % As Xs1 Ys1 are better at minimazation for the middle points 
    % And Xs2 Ys2 are better for the beggining and end parts 
    % Will take best predictor part from each:
    Xs2 = [Xs2(floor(n/2+1):end,:);Xs2(1:floor(n/2),:)];
    Ys2 = [Ys2(floor(n/2+1):end,:);Ys2(1:floor(n/2),:)];
    
    Xs = [Xs2(1:floor(n/4),:); Xs1(floor(n/4+1):floor(3/4*n),:);Xs2(floor(3/4*n+1):n,:)];
    Ys = [Ys2(1:floor(n/4),:); Ys1(floor(n/4+1):floor(3/4*n),:);Ys2(floor(3/4*n+1):n,:)];

%     figure; imshow(orgImg,[]); hold on; plot(Ys, Xs, '--*');
%     figure; imshow(img,[]); hold on; plot(Ys, Xs, '--*');
%     close all;

end
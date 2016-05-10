function startIterations(img, alphas, betas, gammas, templateSize, hoodSize)
    
    % Loading information from template 
    % Starting snake:
    load('landmarks');
    % Neighbors mat:
    load('adjMat');
    % Neighbors distance (1st derivative):
    load('adjDist');
    % Curvature (2nd derivative):
    load('adjCurve');
    
    Xs = landmarks(:,1);
    Ys = landmarks(:,2);
    
    % How many points (landmarks) we have:
    nPoints = size(Xs, 1);
   
    % 
    for p = 1:nPoints
        corrImages(:,:,p) = templateMatching(img, Xs(p), Ys(p), templateSize);
    end
    
    
    for iter=1:3000 %iterations
        
        % Calculating distance between point to next and the average distance:
        snakePointsDist = zeros(nPoints, 1);
        % For each snake point
        for i= 1:nPoints
            % Calculating distance between point and next snake point
            snakePointsDist(i) = sqrt((Xs(i)-Xs(neighborsMat(i,1)))^2 + (Ys(i)-Ys(neighborsMat(i,1)))^2);
        end
        % Computing the average distance from point to point in the current snake
        averageDist = sum(snakePointsDist)/nPoints;
        
        %% For each point in the snake do..
        % This for loop will be changed to randomly pick snake points
        for p=1:nPoints
            xmin = Xs(p);
            ymin = Ys(p);
            
            [normedHoodCont, normedHoodCur, normedHoodTempMatch] = normHood(Xs, Ys, p, hoodSize, averageDist, squeeze(corrImages(:,:,p)));
            
            emin = alphas(p) * normedHoodCont(hoodSize+1, hoodSize+1);
            emin = emin + (beta(p) * normedHoodCur(hoodSize+1, hoodSize+1));
            emin = emin + (gamma(p) * normedHoodEdge(hoodSize+1, hoodSize+1));

            n = hoodSize*2+1;
            newLocation = 0;
            for j = 1:n
                for k = 1:n
                    tempEMin = alphas(p) * normedHoodCont(j,k);
                    tempEMin = tempEMin + (betas(p) * normedHoodCur(j,k));
                    tempEMin = tempEMin + (gammas(p) * normedHoodEdge(j,k));

                    if (tempEMin < emin)
                        emin = tempEMin;
                        newLocation = [j,k];
                    end
                end
            end

            if (newLocation~=0)
                Xs(p) = Xs(p) - ceil(n/2) + newLocation(1);
                Ys(p) = Ys(p) - ceil(n/2) + newLocation(2);
            end
        end  
        
    end
end
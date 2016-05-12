function startIterations(img, alphas, betas, gammas, templateSize, hoodSize, sigma, snakeGif)
    
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
   
    % find the correlation image for each landmark point:
    for p = 1:nPoints
        corrImages(:,:,p) = templateMatching(img, Xs(p), Ys(p), templateSize, sigma);
    end
    
    for iter= 1:3000 %iterations
        if (mod(iter, 100)==0)
            text(50, 50,['iterations=', num2str(iter)],'FontSize',18,'BackgroundColor','black','Color','white');
            scatter(Xs, Ys);
            drawInGif(snakeGif,2);
            figureChildren = get(gca, 'children');
            delete(figureChildren(1:2)); 
        end
        
        % For each point in the snake do..
        % This for loop will be changed to randomly pick snake points
        for p=1:nPoints
            xmin = Xs(p);
            ymin = Ys(p);

            % To get contour (1st derivative) and curvature (2nd '') - Going over all hood points:
            hoodCont = zeros(hoodSize*2+1, hoodSize*2+1);
            hoodCurv = zeros(hoodSize*2+1, hoodSize*2+1);
            countX = 1;
            for x = xmin-hoodSize:xmin+hoodSize
                countY = 1;
                for y = ymin-hoodSize:ymin+hoodSize
                     
                    % Finding the contour and curvature values for the current hood point
                    hoodCont(countX, countY) = calcContEnergy(x, y, p, Xs, Ys, adjMat, adjDist);
                    hoodCurv(countX, countY) = calcCurvEnergy(x, y, p, Xs, Ys, adjMat, adjCurve);

                    countY = countY+1;
                end
                countX = countX+1;
            end
            hoodCorr = corrImages(xmin-hoodSize:xmin+hoodSize, ymin-hoodSize:ymin+hoodSize, p);
            
            hoodCont = hoodCont/max(hoodCont(:));
            hoodCurv = hoodCurv/max(hoodCurv(:));
            
            if (max(hoodCorr(:) ~= 0))
                hoodCorr = -hoodCorr/max(hoodCorr(:));
            end


            emin = alphas * hoodCont(hoodSize+1, hoodSize+1);
            emin = emin + (betas * hoodCurv(hoodSize+1, hoodSize+1));
            emin = emin + (gammas * hoodCorr(hoodSize+1, hoodSize+1));

            n = hoodSize*2+1;
            newLocation = 0;
            for j = 1:n
                for k = 1:n
                    tempEMin = alphas * hoodCont(j,k);
                    tempEMin = tempEMin + (betas * hoodCurv(j,k));
                    tempEMin = tempEMin + (gammas * hoodCorr(j,k));

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
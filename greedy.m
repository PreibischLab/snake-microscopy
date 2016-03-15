%% Function - Greedy snake implementation
% Called by mainSnake

% Based on: Williams & Shah (1992)	

function snake = greedy(img, snake, iterations, hoodSize, orgImg)

    figure;
    
    for iter=1:iterations
        
        imshow(uint8(img));
        impixelinfo;
        hold on;
        title([num2str(iter) , ' iterations, alpha=', num2str(snake(3,1)), ' beta=', num2str(snake(4,1)), ' gamma=', num2str(snake(5,1))])
        plot(snake(1,:), snake(2,:),'LineWidth', 2);
        
        %% Calculating the average snake point to next point distance:
        %  We need this to calculate the contour energy at each point
        
        snakePointsDist = zeros(size(snake,2),1);
        % For each snake point
        for i= 1:(size(snake,2)-1)
            % Calculating distance between point and next snake point
            snakePointsDist(i) = sqrt((snake(1,i)-snake(1,i+1))^2 + (snake(2,i)-snake(2,i+1))^2);
        end
        % Calculating distance between last and first point in snake
        i = i+1;
        snakePointsDist(i) = sqrt((snake(1,i)-snake(1,1))^2 + (snake(2,i)-snake(2,1))^2);

        % Computing the average distance from point to point in the current snake
        averageDist = sum(snakePointsDist)/size(snake,2);

        %% For each point in the snake do..
        % If this for loop will be changed to randomly picking snake points
        % there is less of a chance the snake will "spin"
        for i=1:size(snake,2)

            xmin = snake(1,i);
            ymin = snake(2,i);

            [normedHoodCont, normedHoodCur, normedHoodEdge] = normHood(img, snake, i, hoodSize, averageDist);

            emin = snake(3,i) * normedHoodCont(hoodSize+1, hoodSize+1);
            emin = emin + (snake(4,i) * normedHoodCur(hoodSize+1, hoodSize+1));
            emin = emin + (snake(5,i) * normedHoodEdge(hoodSize+1, hoodSize+1));

            n = hoodSize*2+1;
            newLocation = 0;
            for j = 1:n
                for k = 1:n
                    tempEMin = snake(3,i) * normedHoodCont(j,k);
                    tempEMin = tempEMin + (snake(4,i) * normedHoodCur(j,k));
                    tempEMin = tempEMin + (snake(5,i) * normedHoodEdge(j,k));

                    if (tempEMin < emin)
                        emin = tempEMin;
                        newLocation = [j,k];
                    end
                end
            end

            if (newLocation~=0)
                snake(1,i) = snake(1,i) - ceil(n/2) + newLocation(1);
                snake(2,i) = snake(2,i) - ceil(n/2) + newLocation(2);
            end

        end
        
    end
end
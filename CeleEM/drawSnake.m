function drawSnake()
    
    img = imread('../microscopyImages/z=12.png');
    img = imresize(img,0.4);
    
    [img1, rect] = imcrop(img);
    clear('img');
    pause();
    
    figure; imshow(img1); hold on;
    
    snakeGif = 'snakeGif2.gif';
    drawInGif(snakeGif, 1);
    
    for i = 2:2
        
        load([num2str(i) '.mat']);
        for j= 1:30
            text(20, 20,'Balloon Snake','FontSize',8,'BackgroundColor','black','Color','white');
            text(20, 50,['iterations=', num2str(j*100)],'FontSize',8,'BackgroundColor','black','Color','white');
            plot(cell2mat(Xss(j))-rect(1), cell2mat(Yss(j))-rect(2), 'm', 'LineWidth', 2);
            drawInGif(snakeGif, 2);
            figureChildren = get(gca, 'children');
            delete(figureChildren(1:2)); 
        end
        plot(cell2mat(Xss(j)), cell2mat(Yss(j)), 'm', 'LineWidth', 2);
    end
    
    for i = 2:2
        
        load(['1' num2str(i) '.mat']);
        s = size(Xss,1);
        for j= 1:30
            text(20, 20,'Shrinking Snake','FontSize',8,'BackgroundColor','black','Color','white');
            text(20, 50,['iterations=', num2str(j*100)],'FontSize',8,'BackgroundColor','black','Color','white');
            plot(cell2mat(Xss(j))-rect(1), cell2mat(Yss(j))-rect(2), 'b', 'LineWidth', 2);
            drawInGif(snakeGif,2);
            figureChildren = get(gca, 'children');
            delete(figureChildren(1:2)); 
        end
        plot(cell2mat(Xss(j)), cell2mat(Yss(j)), 'b', 'LineWidth', 2);
    end
    
    figureChildren = get(gca, 'children');
    delete(figureChildren(1:end-1)); 
    for i = 2:2
        
        load(['3' num2str(i) '.mat']);

        text(20, 20,'Viterbi Algorithm','FontSize',8,'BackgroundColor','black','Color','white');
        plot(Xs-rect(1), Ys-rect(2), 'g', 'LineWidth', 2);
        drawInGif(snakeGif,2);
    end


end
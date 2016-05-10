function [Xs, Ys] = snakeIterations(img, Xs, Ys, alphas, betas, delta, rho, snakeGif)
    
    % Find the X and Y derivatives of the img (img already after gaussian,sobel,normalization)
    edgeImgX = imfilter(img, fspecial('sobel')');
    edgeImgY = imfilter(img, fspecial('sobel'));
    
    % Build coefficient pentadiagonal matrix - alpha and beta 
    A = buildCoMatA(alphas, betas, delta, size(Xs,1));
    
    if (rho == 0)
        NsX = 0;
        NsY = 0;
    end
    
    Xss = cell(50,1);
    Yss = cell(50,1);
    
    %% Start iterations on snake
    for i = 1:6000%iterations
        
        % Making a gif
        if (mod(i,100) == 0)
            
            Xss(i/100,1) = {Xs};
            Yss(i/100,1) = {Ys};
            
%             text(50, 50,['iterations=', num2str(i)],'FontSize',18,'BackgroundColor','black','Color','white');
%             plot(Xs, Ys, 'm', 'LineWidth', 2);
%             drawInGif(snakeGif,2);
%             figureChildren = get(gca, 'children');
%             delete(figureChildren(1:2)); 
        end

        % Every few iterations to resample the snake points to keep equal distance
        if (mod(i,100) == 0)
            
            [Xs, Ys] = resampleSnake(img, Xs, Ys);

            % Build A again - because n might change
            A = buildCoMatA(alphas, betas, delta, size(Xs,1));
           
        end
        
        % For balloon snake:
        if (rho ~= 0)
            [NsX, NsY] = addExpandingForce(Xs,Ys);
        end
        
        % Find image values at snake points with interpolate 
        imgX = interp2(double(edgeImgX), Xs, Ys, '*linear') + rho * NsX;
        imgY = interp2(double(edgeImgY), Xs, Ys, '*linear') + rho * NsY;
        
        % New snake:
        Xs = A \ ((1/delta)*Xs + imgX);
        Ys = A \ ((1/delta)*Ys + imgY);
        
        if (any(isnan(Xs)))
            error('Out Of Bounds!! :( ');
        end
        
    end
    
    save('40.mat','Xss','Yss');
    plot(Xs, Ys, 'm', 'LineWidth', 2);
    
end
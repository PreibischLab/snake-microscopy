function drawInGif(snakeGif, isFirst)    
    drawnow;
    im = frame2im(getframe(1));
    [imind, cm] = rgb2ind(im,256);
    if (isFirst==1)
        imwrite(imind, cm, snakeGif,'gif', 'Loopcount', inf);
    else
        imwrite(imind, cm, snakeGif, 'gif', 'WriteMode', 'append');
    end
end
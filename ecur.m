%% Evaluating Contour Curvatur - "The second order differential"
function curve = ecur(x, y, s, snake)

    if ((s>1) && (s<size(snake,2)))
        curve = (snake(1,s-1) - 2*x + snake(1,s+1))^2 + (snake(2,s-1) - 2*y + snake(2,s+1))^2;
    elseif (s==1)
        curve = (snake(1,size(snake,2)) - 2*x + snake(1,2))^2 + (snake(2,size(snake,2)) - 2*y + snake(2,2))^2;
    else
        curve = (snake(1,size(snake,2)-1) - 2*x + snake(1,1))^2 + (snake(2,size(snake,2)-1) - 2*y + snake(2,1))^2;
    end
end
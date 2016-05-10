function A = buildCoMatA(alpha, beta, delta, n)

    % Defining the spacing between points in snake
    % Probably irrelevant - hence h=1
    h = 1;
    
    % Get number of snake points:
    alphas(1:n) = alpha;
    betas(1:n) = beta;
    
    % Since the last point of the snake is the point before the first and so on
    % Defining the matching alphas and betas to those cases at the edges:
    alphasPrev = [alphas(n) alphas(1:n-1)];
    betasPrev = [betas(n) betas(1:n-1)];
    betasNext = [betas(2:n) betas(1)];
    
    % Those are just the discrete approximation of the snake equations:
    a = betasNext/h^4;
    b = -2*(betas+betasNext)/h^4 - alphas/h^2;
    c = (betasPrev+4*betas+betasNext)/h^4 + (alphasPrev+alphas)/h^2;
    d = -2*(betasPrev+betas)/h^4 - alphasPrev/h^2;
    e = betasPrev/h^4;
    
    % Building the matrix from those equations:
    
    % Setting the values we want in each non-zero diagonal:
    temp = [d; e; a; b; c; d; e; a; b]';
    A = full(spdiags(temp, [-(n-1) -(n-2) -2 -1 0 1 2 (n-2) (n-1)], n, n));
    
    % Detla controls step size:
    A = A + (1/delta)*eye(n);
end
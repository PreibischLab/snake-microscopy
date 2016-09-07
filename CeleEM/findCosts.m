function E = findCosts(X, Y, Eext, i, lambda)

    m = size(X, 2);
    
    % Init arrays to hold internal cost of step
    Eint = nan(m, m, m);
    E = nan(m, m, m);
    
    for prev = 1:m
        for curr = 1:m
            for next = 1:m

                nomX = X(i+1, prev) - 2*X(i, curr) + X(i-1, next);
                nomY = Y(i+1, prev) - 2*Y(i, curr) + Y(i-1, next);

                denX = X(i+1, next) - X(i-1, prev);
                denY = Y(i+1, next) - Y(i-1, prev);

                nom = sqrt(nomX.^2 + nomY.^2);
                den = sqrt(denX.^2 + denY.^2);

                Eint(prev, next, curr) = (nom ./ den).^2;
            end
        end
    end
    
    for curr = 1:m
        E(:,:,curr) = lambda * Eint(:,:,curr) + (1 - lambda) * Eext(curr);
    end
    
end
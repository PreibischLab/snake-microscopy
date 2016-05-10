function [NsX, NsY] = addExpandingForce(Xs, Ys)

    n = size(Xs,1);
    
    % Initialize expanding force X and Y arrays
    DsX = nan(n,1);
    DsY = nan(n,1);

    for i = 1:n
        if (i<n)
            % Calcs the distance between a snake point to the next
            DsX(i) = -(Ys(i+1) - Ys(i));
            DsY(i) = (Xs(i+1) - Xs(i));
        else
            % For the Sn point - S1 is the next point
            DsX(n) = -(Ys(1) - Ys(n));
            DsY(n) = (Xs(1) - Xs(n));
        end
    end
    
    % calcs the actual distance between each 2 points
    magnitude = sqrt((DsX.^2) + (DsY.^2));
    
    if (any(isnan(magnitude)))
        magnitude
        error('magnitue error in finding norm');
    end
    
    % Since we want the norm - divide by the distance
    NsX = - DsX ./ magnitude;
    NsY = - DsY ./ magnitude;
    
end
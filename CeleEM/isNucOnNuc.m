% This script gets a nuc from one slice and another slice:

% It checks to see if there is a nucleus on the area of this nuc in another
% slice - if so - we assume its the same nuc (assuming that the slices are close to one another)

function nucNumInOtherSlice = isNucOnNuc(nuc, sliceNum)
    
    Xs = nuc(:,1);
    Ys = nuc(:,2);
    load('slicesNucs');
    
    nucs = slicesNucs(:,sliceNum);
    nucs = nucs(~cellfun('isempty',nucs));
    
    xStart = Xs(1)-500;
    yStart = Ys(1)-500;
	
    Xcontour = round(Xs(:)-xStart);
	Ycontour = round(Ys(:)-yStart);
    
    labelCurr = double(poly2mask(Ycontour, Xcontour, 1000, 1000));
    
    % If we don't find a nuc there, this variable will stay nan:
    nucNumInOtherSlice = nan;
    
    counter = 1;
    
    for i=1:size(nucs,1)
        if (~isnan(nucs{i}(1)))            
            if ((abs(nucs{i}(1,1) - Xs(1)) < 150) && (abs(nucs{i}(1,2) - Ys(1)) < 300))

                contour = nucs{i};

                Xcontour = round(contour(:,1)-xStart);
                Ycontour = round(contour(:,2)-yStart);
                labelNext = double(poly2mask(Ycontour, Xcontour, 1000, 1000));

                isSameNuc = max(max(labelCurr+labelNext));
                if (isSameNuc==2)
                    % The function returns the current number of the nuc in the slice
                    nucNumInOtherSlice(counter) = i;
                    counter = counter+1;
                end
            end
        end
    end

end
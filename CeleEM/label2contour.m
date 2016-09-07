% Label to contour gets a mask image (annotated image) of all nuclei in slice and turns it into
% a cell array - each cell is one nuclous, defined by its contour (Xs and Ys)

function nucsContours = label2contour(labelImg, sample)
    
    sizeLabel = size(labelImg);
    
    % Finds all the different objects (components) in the image
    % So for each nucleous we get the location of all of its pixels.
    nucs = bwconncomp(labelImg);
    
    % Get count how many nucleous in this slice:
    nNucs = nucs.NumObjects;
    % Initialize cell array to keep all nuclei contours
    nucsContours = cell(nNucs,1);
    
    for i=1:nNucs
        % Assign to nuc all of its pixels
        nuc = nucs.PixelIdxList{i};
        % Initialize an image to put nuc in
        polyg = zeros(sizeLabel);
        % Make mask with only nuc in
        % This is only possible in the linear index and not in col and row sub, so leaving nuc untouched works.
        polyg(nuc) = 1;

        % Find first contour pixel in nuc
        [tempX, tempY] = find(polyg,1);
        % Trace the contour
        edgePolyg = bwtraceboundary(polyg, [tempX tempY],'N');
        % Delete the last pixel in contour (since it's equal to first)
        edgePolyg = edgePolyg(1:end-1,:);
        % To mark contour only by one in x (defined by sample) pixels:
%         nPoints = size(edgePolyg,1) / 5;
%         sample = floor(size(edgePolyg,1)/nPoints);
        nucContour = [edgePolyg(1:sample:end, 1), edgePolyg(1:sample:end, 2)];
        nucContour(end,:) = [];

        % Assign nucleous contour values
        nucsContours{i} = nucContour;
    end

end
% [ldmkMoving, out] = startIterations(img, shapeModel, InitialLandmarks, templates, hoodSize, ldmkIdx)
% Parameters:
%   img:
%     - image on which the landmarks should be retrived
% 
%   shapeModel:
%     - structure containing :
%        - paramPerLdmk: nLdmk x 3 matrix
%             row: [alphas betas gammas]
%               Alpha - weight of elasticity for list of points (distance)
%               Beta  - weight of curvature (abs angular difference)
%               Gamma - weight of ext force
%        - adjIndices: nLdmk x 3 matrix 
%              index of adjascent landmarks
%        - adjCurve: nLdmk x 2 matrix 
%              oriented angle between L and Adj1 and Adj2
%                                     L and Adj2 and Adj3
%        - adjDist: nLdmk x 3 matrix 
%              distance between landmarks and their adjascent points
% 
%   InitialLandmarks: 
%        - initial landmark position
%           if no better, use the position of landmarks on the template
% 
%   templates:
%         - k x k x nLdmk matrix 
%              stack images  of the template neighborhood that are
%              matched through cross correlation 
% 
%   hoodSize:  
%         - size of the neihborhood on which the calculation are made
%                at each step
% 
%   ldmkIdx (optional): 
%         - index of landmark that will be updated by the snake.
%               By default, all landmarks are updated.
%                



function [ldmkMoving, out, ldmkHistory] = startIterations(img, shapeModel, templates, hoodSize, ldmkIdx, userStartPosition)

boolOutputVis = nargout > 1;
% strVis = 'last';
strVis = 'init';

plotEachLandmark = false(1);

%% "hidden" parameters
nIterMax = 10;
sigma = 3;
bigHood = 150;

%% get the shape model information
ldmkClass        = shapeModel.ldmkClass;
segments         = shapeModel.segments;
relPosInSegments = shapeModel.relPosInSegments;
weights          = shapeModel.weights;

% adjacency is used only for curvature
adjIdx          = shapeModel.adjIndices;
adjCurveVal     = shapeModel.adjCurveVal;
adjCurveWeight  = shapeModel.adjCurveWeight;

distIdx         = shapeModel.distIdx;
distSubVal      = shapeModel.distSubVal;
distSubWeight   = shapeModel.distSubWeight;

curvIdx         = shapeModel.curvIdx;
curvSubVal      = shapeModel.curvSubVal;
curvSubWeight   = shapeModel.curvSubWeight;

processOrder    = shapeModel.processOrder;
% centroidAvDist  = shapeModel.centroidAvDist;

if (exist('userStartPosition','var') && ~isempty(userStartPosition))
    startPosition = round(userStartPosition);
else
    startPosition = round(shapeModel.startPosition);
end

% normalize the relative weight
adjCurveWeight = cellfun(@(x) x/sum(x(:)),adjCurveWeight,'Unif',false);
distSubWeight  = cellfun(@(x) x/sum(x(:)),distSubWeight ,'Unif',false);
curvSubWeight  = cellfun(@(x) x/sum(x(:)),curvSubWeight ,'Unif',false);

%% checks

if any(cellfun(@(x) sum(x<=0 & x>=1)>0, relPosInSegments))
    error('relative distance must be in the interval ]0 1[')
end

idxSemi    = reshape(find(ldmkClass==3),1,[]);
ixSegment  = zeros(size(startPosition,1),1);
segmentAdj = zeros(size(startPosition,1),2); 
segmentAdjrelDist = zeros(size(startPosition,1),2); 
for i=idxSemi
    % ind its segment
    A = cellfun(@(Set) ismember(i,Set), segments);
    if sum(A)~=1
        error('semi attributed to more than one segment');
    end
    ixSegment(i) = find(A);
    
    % position within the segment
    j = find(segments{ixSegment(i)}==i);
    
    % neighbour
    segmentAdj(i,:) = [ segments{ixSegment(i)}(j-1),...
                        segments{ixSegment(i)}(j+1) ];
    R = relPosInSegments{ixSegment(i)};
    segmentAdjrelDist(i,:) = [R(j) - R(j-1), R(j+1) - R(j)]; 
end

%% deals with the inputs
% get the number of landmaks
nLdmk = size(startPosition, 1);

startPosition = round(startPosition);
ldmkMoving = startPosition;

if (~exist('ldmkIdx','var')) || isempty(ldmkIdx)
    ldmkIdx = 1:nLdmk;
end

% reshape for the for loop
ldmkIdx = reshape(ldmkIdx,1,[]);

[i,locb] = ismember(processOrder,ldmkIdx);
idxLdmkOrder = reshape(ldmkIdx(locb(i)),1,[]);

% add non processed landmarks at the end
idxLdmkOrderFull = [idxLdmkOrder setdiff(1:nLdmk,idxLdmkOrder)];

%%  precalculate the correlation image for each landmark point:
% big neighborhood of each point: max searching area (for sake of speed)
xIdx = zeros(nLdmk, bigHood*2+1);
ixIdx = zeros(nLdmk, bigHood*2+1);
yIdx = zeros(nLdmk, bigHood*2+1);
iyIdx = zeros(nLdmk, bigHood*2+1);
for p = ldmkIdx
    xIdx(p,:) = ((-bigHood):(bigHood)) + ldmkMoving(p,1);
    yIdx(p,:) = ((-bigHood):(bigHood)) + ldmkMoving(p,2);
end
sbigHood = size(xIdx,2);

% calculate the cross correlation on the smaller windows
corrImages = -ones(size(yIdx,2),size(xIdx,2),nLdmk);
G = fspecial('gaussian',[5 5],sigma);
imgFilter = imfilter(img,G,'same');
for p = 1:nLdmk
    iX = (xIdx(p,:)>0) & (xIdx(p,:)<=size(img,1));
    ixIdx(p,:) = iX;
    iY = (yIdx(p,:)>0) & (yIdx(p,:)<=size(img,2));
    iyIdx(p,:) = iY;
    corrImages(iY,iX,p) = normxcorr2e(templates{p}, imgFilter(yIdx(p,iX),xIdx(p,iY)), 'same');
end

% normalize between 0 and 1 (0 best match, 1 worse match)
corrImages = 1-(corrImages+1)/2;

%% variable initialization and precalculations
% initialize neighborhood matrices
hoodDist = zeros(hoodSize*2+1, hoodSize*2+1);
hoodCurv = zeros(hoodSize*2+1, hoodSize*2+1);
hoodCorr = zeros(hoodSize*2+1, hoodSize*2+1);

% precompute the coordinates of neighborhood pixels
[X, Y] = meshgrid((-hoodSize):(hoodSize), (-hoodSize):(hoodSize));
hoodCoord = [X(:) Y(:)];

% precompute the distance to center used in case of end points landmarks
hoodDistPre  =   abs(complex(X,Y));

%  hoodAngle = angle(complex(X,Y));

% gaussian filter for the
hg = fspecial('gaussian', [3 3], 2);
%%
if boolOutputVis
    out.metric    = nan(size(img));
    out.distance  = nan(size(img));
    out.curvature = nan(size(img));
    out.total     = nan(size(img));
else
    out.metric    = [];
    out.distance  = [];
    out.curvature = [];
    out.total     = [];
end

%% iteration loop
ldmkHistory = cell(nIterMax+1,1);
ldmkHistory{1} = ldmkMoving;
for iter= 1:nIterMax %iterations
    if boolOutputVis && strcmp(strVis,'last')
        out.metric{1}    = nan(size(img));
        out.distance{1}  = nan(size(img));
        out.curvature{1} = nan(size(img));
        out.total{1}     = nan(size(img));
    end
    
    % estimate the centroid size: average distance to the center
    centroidAvDist = mean(sqrt(sum(bsxfun(@minus, ldmkMoving, mean(ldmkMoving,1)).^2,2)));
    
    ldmkMovingOld = ldmkMoving;
    for p = idxLdmkOrderFull
        
        HCoordTmp = bsxfun(@plus, hoodCoord, ldmkMoving(p,:));
        
        %% distance
        % pair distance between adjascent points and all neighboring points
        if ldmkClass(p) == 1 || ldmkClass(p) == 2 % for cross points and end points
            
            adjLdmkTmp = ldmkMoving(distIdx{p},:);
            
             % get the  distance relative to centroid size
            D = pdist2(HCoordTmp, adjLdmkTmp) / centroidAvDist;
            
            % difference to the expected one
            D = bsxfun(@minus, D , distSubVal{p});
            
            % average the distance
            hoodDist(:) = mean(bsxfun(@times, abs(D), distSubWeight{p}),2);
            
        else % for semi landmarks
            
            % overall distance of the segment
            Dtmp = squareform(pdist(ldmkMoving(segments{ixSegment(p)},:)));
            sgmDist = sum(Dtmp(diag(ones(1,numel(segments{ixSegment(p)})-1),1)>0));
            
            % distance to neighbor
            adjLdmkTmp = ldmkMoving(segmentAdj(p,:),:);
            D = pdist2(HCoordTmp, adjLdmkTmp);
            D = bsxfun(@minus, D , segmentAdjrelDist(p,:)*sgmDist );
            hoodDist(:) = mean(abs(D(:,1)),2);
        end
        
        %% curvature
        if ldmkClass(p) ==1 || ldmkClass(p) ==2
            hoodCurv(:) = 0;
            if ~isempty(curvIdx{p})
                for n = 1:size(curvIdx{p},1)
                    adjLdmkTmp = permute(ldmkMoving(curvIdx{p}(n,:),:), [3 2 1]);
                    hoodCurv(:) = hoodCurv(:) + ...
                        f_AgDiff(HCoordTmp, adjLdmkTmp, curvSubVal{p}(n)) * curvSubWeight{p}(n);
                end
            else
                hoodCurv(:) = hoodDistPre(:);
            end
            
        else
            % calculation for semi ldmk only based on adjacency neighborhood
            hoodCurv(:) = 0;
            for n = 1:size(adjIdx{p},1)
                adjLdmkTmp = permute(ldmkMoving(adjIdx{p}(n,:),:), [3 2 1]);
                hoodCurv(:) = hoodCurv(:) + ...
                    f_AgDiff(HCoordTmp, adjLdmkTmp, adjCurveVal{p}(n)) * adjCurveWeight{p}(n);
            end
            
        end
        
        %% find minimum
        % normalize curvature and distance
        hoodDist = hoodDist/max(hoodDist(:));
        hoodCurv = hoodCurv/max(hoodCurv(:));

        
        % coordinates of the square windows around the landmark
        xWindows = X(1,:) + ldmkMoving(p,1) - xIdx(p,1);
        yWindows = Y(:,1) + ldmkMoving(p,2) - yIdx(p,1);
        iX = (xWindows > 0) & (xWindows <= sbigHood);
        iY = (yWindows > 0) & (yWindows <= sbigHood);
        
        hoodCorr(:) = 1;
        % do not normalize the correlation because bad match should not
        % drive the point movment too much
        hoodCorr(iY,iX) = corrImages(yWindows(iY), xWindows(iX), p);

        % final matrix
        tmp = (weights(p,1) * hoodDist + ...
               weights(p,2) * hoodCurv + ...
               weights(p,3) * hoodCorr)./sum(weights(p,:));

        tmp = filter2(hg,tmp, 'same');
        m = max(tmp(:));
        tmp([1 end],1:end) = m;
        tmp(1:end,[1 end]) = m;

        
        if boolOutputVis && ((strcmp(strVis,'init') && iter ==1) || strcmp(strVis,'last'))
            Xtmp = X(1,:) + ldmkMoving(p,1);
            Ytmp = Y(:,1) + ldmkMoving(p,2);
            
            out.metric(Ytmp,Xtmp)    = hoodCorr;
            out.distance(Ytmp,Xtmp)  = hoodDist;
            out.curvature(Ytmp,Xtmp) = hoodCurv;
            out.total(Ytmp,Xtmp)     = tmp;
        end

        % find the best position and move the landmark
        if ismember(p,idxLdmkOrder)
            [~,iMin] = min(tmp(:));
            ldmkMoving(p,:) = ldmkMoving(p,:) + [X(iMin), Y(iMin)];
        end
        
        %% plot (optional)
        if plotEachLandmark %&& (p==1) % || p ==14)
            xWindows = X(1,:) + ldmkMovingOld(p,1);
            yWindows = Y(:,1) + ldmkMovingOld(p,2);
            current = ldmkMovingOld(p,:);
            
            subplot(2,3,1)
            imagesc(xWindows, yWindows, tmp)
            hold on
            scatter(current(1), current(2),'+','w');
            scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+','r');
            hold off
            axis equal tight
            title('total')

            subplot(2,3,2)
            imagesc(xWindows, yWindows, hoodCurv)
            hold on
            scatter(current(1), current(2),'+','w');
            scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+','r');
            hold off
            axis equal tight
            title('curv')

            subplot(2,3,3)
            imagesc(xWindows, yWindows, hoodDist)
            hold on
            scatter(current(1), current(2),'+','w');
            scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+','r');
            hold off
            axis equal tight
            title('dist')

            subplot(2,3,4)
            imagesc(xWindows, yWindows, hoodCorr)
            hold on
            scatter(current(1), current(2),'+','w');
            scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+','r');
            hold off
            axis equal tight
            title('corr')

            subplot(2,3,5)
            imagesc(xWindows, yWindows, imgFilter(yWindows,xWindows))
            hold on
            scatter(current(1), current(2),'+','g');
            scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+','r');
            hold off
            axis equal tight
            title('image')

            subplot(2,3,6)
            imagesc(templates{p})
            hold on
            scatter(size(templates{p},2)/2, size(templates{p},1)/2,'+','r');
            hold off
            axis equal tight
            title('template')
        end
    end 
    ldmkHistory{iter+1} = ldmkMoving;

    % stop criterion
    if sum(sqrt(sum((ldmkMovingOld - ldmkMoving ).^2,2))) < 1
        break
    end
    
end

ldmkHistory = cell2mat(reshape(ldmkHistory(1:iter+1),1,1,[]));
    
    
    
end

function anglDiff = f_AgDiff(HoodCoordinate, adjascentPos, expectedAngle)

% vector between ldmk and adjascent ldmks
V = bsxfun(@plus, -HoodCoordinate, adjascentPos);

% calulate the oriented angle bt these
% dot and cross product
vDot = sum(prod(V,3),2);
vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
    
% get the angle
anglDiff = atan2(vCross, vDot);
    
% warp it to that the results are comparible
anglDiff = fAngularDiff(anglDiff, expectedAngle);

end






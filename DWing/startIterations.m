% ldmkMoving = startIterations(img, shapeModel, InitialLandmarks, templates, hoodSize, ldmkIdx)
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



function ldmkMoving = startIterations(img, shapeModel, InitialLandmarks, templates, hoodSize, ldmkIdx)

%% "hidden" parameters
sigma = 3;
bigHood = 100;
plotEachLandmark = false(1);

%% get the shape model information
paramPerLdmk = shapeModel.paramPerLdmk;
adjIndices   = shapeModel.adjIndices;
adjCurve     = shapeModel.adjCurvature;
adjDist      = shapeModel.adjDistance;

%% deals with the inputs
% get the number of landmaks
nLdmk = size(InitialLandmarks, 1);

InitialLandmarks = round(InitialLandmarks);
ldmkMoving = InitialLandmarks;

if (~exist('ldmkIdx','var')) || isempty(ldmkIdx)
    ldmkIdx = 1:nLdmk;
end

% reshape for the for loop
ldmkIdx = reshape(ldmkIdx,1,[]);


% refer to NaN index in adjascent point as index to landmark of NaN
% coordinates (facilitate calculations)
ldmkMoving = [ldmkMoving; NaN NaN];
adjIndices(isnan(adjIndices)) = nLdmk+1;


%%  precalculate the correlation image for each landmark point:
% big neighborhood of each point: max searching area (for sake of speed)
xIdx = zeros(nLdmk, bigHood*2+1);
yIdx = zeros(nLdmk, bigHood*2+1);
for p=1:nLdmk
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
    iY = (yIdx(p,:)>0) & (yIdx(p,:)<=size(img,2));
    corrImages(iY,iX,p) = normxcorr2e(templates(:,:,p), imgFilter(yIdx(p,iX),xIdx(p,iY)), 'same');
end

% normalize between 0 and 1 (0 best match, 1 worse match)
corrImages = 1-(corrImages+1)/2;

%% variable initialization and precalculations
% initialize neighborhood matrices
hoodCont = zeros(hoodSize*2+1, hoodSize*2+1);
hoodCurv = zeros(hoodSize*2+1, hoodSize*2+1);
hoodCorr = zeros(hoodSize*2+1, hoodSize*2+1);

% precompute the coordinates of neighborhood pixels
[X, Y] = meshgrid((-hoodSize):(hoodSize), (-hoodSize):(hoodSize));
hoodCoord = [X(:) Y(:)];
NNeighbor = sum(adjIndices~=nLdmk+1,2);

% precompute the distance to center used in case of end points landmarks
hoodDist  =   abs(complex(X,Y));

% hoodAngle = angle(complex(X,Y));

% gaussian filter for the
hg = fspecial('gaussian', [3 3], 2);

%% iteration loop
for iter= 1:20 %iterations
    ldmkMovingOld = ldmkMoving;
    for p = ldmkIdx
        if iter==2
            1;
        end
        hoodCoordTmp = bsxfun(@plus, hoodCoord, ldmkMoving(p,:));
        adjLdmkTmp = ldmkMoving(adjIndices(p,:),:);
        
        % pair distance between adjascent points and all neighboring points
        D = pdist2(hoodCoordTmp, adjLdmkTmp);
        D = bsxfun(@minus, D , adjDist(p,:));
        D(isnan(D)) = 0;
        
        % average the distance
        hoodCont(:) = sum(abs(D),2)/NNeighbor(p);

        % different calculations depending on the number of neighbor
        adjLdmkTmpPerm= permute(adjLdmkTmp(:,:),[3 2 1]);
        switch NNeighbor(p)
            case 1
                % to improve ? perpendicular to the vector ?
                hoodCurv(:) = hoodDist(:);

            case 2
                % vector between ldmk and adjascent ldmks
                V = bsxfun(@plus, -hoodCoordTmp, adjLdmkTmpPerm(1,:,[1 2]));
            % calulate the oriented angle bt these
                % dot and corss product
                vDot = sum(prod(V,3),2);
                vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
                % get the angle
                tmp = atan2(vCross, vDot);

                % warp it to that the results are comparible
                tmp = mod(tmp- adjCurve(p,1),2*pi);
                hoodCurv(:) = min(2*pi-tmp,tmp);

            case 3
                V = bsxfun(@plus, -hoodCoordTmp, adjLdmkTmpPerm(1,:,[1 2]));
                vDot = sum(prod(V,3),2);
                vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
                tmp = atan2(vCross, vDot);
                tmp = mod(tmp- adjCurve(p,1),2*pi);

                V = bsxfun(@plus, -hoodCoordTmp, adjLdmkTmpPerm(1,:,[2 3]));
                vDot = sum(prod(V,3),2);
                vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
                tmp2 = atan2(vCross, vDot);
                tmp2 = mod(tmp2- adjCurve(p,2),2*pi);

                % averagea of the two angles
                hoodCurv(:) = (min(2*pi-tmp ,tmp ) + ...
                               min(2*pi-tmp2,tmp2)   )/2;
            otherwise
                error('oups')
        end

        % normalize curvature and distance
        hoodCont = hoodCont/max(hoodCont(:));
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
        tmp = paramPerLdmk(p,1) * hoodCont + ...
              paramPerLdmk(p,2) * hoodCurv + ...
              paramPerLdmk(p,3) * hoodCorr;

        tmp = filter2(hg,tmp, 'same');
        m = max(tmp(:));
        tmp([1 end],1:end) = m;
        tmp(1:end,[1 end]) = m;

        % find the best position and move the landmark
        [~,iMin] = min(tmp(:));
        ldmkMoving(p,:) = ldmkMoving(p,:) + [X(iMin), Y(iMin)];
        
        %% plot (optional)
        if plotEachLandmark%  && (p==12) % || p ==14)
            xWindows = X(1,:) + ldmkMovingOld(p,1);
            yWindows = Y(:,1) + ldmkMovingOld(p,2);
            current = ldmkMovingOld(p,:);
            figure(2)
            imagesc(xWindows, yWindows, tmp)
            hold on
            scatter(current(1), current(2),'+','w');
            scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+','r');
            hold off
            axis equal tight
            title('total')

            figure(3)
            imagesc(xWindows, yWindows, hoodCurv)
            hold on
            scatter(current(1), current(2),'+','w');
            scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+','r');
            hold off
            axis equal tight
            title('curv')

            figure(4)
            imagesc(xWindows, yWindows, hoodCont)
            hold on
            scatter(current(1), current(2),'+','w');
            scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+','r');
            hold off
            axis equal tight
            title('dist')

            figure(5)
            imagesc(xWindows, yWindows, hoodCorr)
            hold on
            scatter(current(1), current(2),'+','w');
            scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+','r');
            hold off
            axis equal tight
            title('corr')

            figure(6)
            imagesc(xWindows, yWindows, imgFilter(yWindows,xWindows))
            hold on
            scatter(current(1), current(2),'+','g');
            scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+','r');
            hold off
            axis equal tight
            title('image')

            figure(7)
            imagesc(templates(:,:,p))
            hold on
            scatter(16, 16,'+','r');
            hold off
            axis equal tight
            title('template')
        end
    end 


    % stop criterion
    
    if sum(sqrt(sum((ldmkMovingOld(1:end-1,:)-...
                     ldmkMoving(1:end-1,:)  ).^2,2))) < 1
        break
    end
end
    
ldmkMoving = ldmkMoving(1:end-1,:);
    
    
end








function ldmkMoving = startIterations(img, ldmk, adjM, paramPerLdmk, hoodSize, templates, ldmkIdx, ldmkMoving)
sigma = 3;
ldmk = round(ldmk);
if (~exist('ldmkIdx','var')) || isempty(ldmkIdx)
    ldmkIdx = 1:size(ldmk,1);
end
if (~exist('ldmkMoving','var')) || isempty(ldmkMoving)
    ldmkMoving = ldmk;
end
ldmkMoving = round(ldmkMoving);

% reshape for the for loop
ldmkIdx = reshape(ldmkIdx,1,[]);

nLdmk = size(ldmk, 1);
ldmk = [ldmk; NaN NaN];
ldmkMoving = [ldmkMoving; NaN NaN];
adjM(isnan(adjM)) = nLdmk+1;

% table with dim1 = ldmk index,
%            dim2 = [X Y] coordinates 
%            dim3 = adj index,
% for easier calculations
adjLdmk = reshape(ldmk(adjM,:),size(adjM,1),size(adjM,2),2);
adjLdmk = permute(adjLdmk,[1 3 2]);

%% estimate distance between adjascent points
% distance to adjascent points
adjDist = sqrt( sum( bsxfun(@minus, adjLdmk, ldmk(1:end-1,:)).^2 ,2) );
adjDist = permute(adjDist,[1 3 2]);

adjCurve = nan(size(adjLdmk,1),2);
% curvature (not working ?)
%     adjCurve(:,1) =  sqrt(sum( (sum(adjLdmk(:,:,[1 2]),3) - 2*ldmk(iLdmk,:)).^2 ,2));
%     adjCurve(:,2) =  sqrt(sum( (sum(adjLdmk(:,:,[1 3]),3) - 2*ldmk(iLdmk,:)).^2 ,2));
%     adjCurve(:,3) =  sqrt(sum( (sum(adjLdmk(:,:,[2 3]),3) - 2*ldmk(iLdmk,:)).^2 ,2));
%     adjCurve(isnan(adjCurve)) = 0;

% oriented angle
% (https://stackoverflow.com/questions/3486172/angle-between-3-points)
V = bsxfun(@minus, adjLdmk(:,:,[1 2]), ldmk(1:end-1,:));
vDot = sum(prod(V,3),2);
vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
adjCurve(:,1) = atan2(vCross, vDot);


V = bsxfun(@minus, adjLdmk(:,:,[2 3]), ldmk(1:end-1,:));
vDot = sum(prod(V,3),2);
vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
adjCurve(:,2) = atan2(vCross, vDot);

%% big neighborhood of each point: max searching area (for sake of speed)
bigHood = 100;
xIdx = zeros(nLdmk, bigHood*2+1);
yIdx = zeros(nLdmk, bigHood*2+1);
for p=1:nLdmk
    xIdx(p,:) = ((-bigHood):(bigHood)) + ldmkMoving(p,1);
    yIdx(p,:) = ((-bigHood):(bigHood)) + ldmkMoving(p,2);
end
sbigHood = size(xIdx,2);
%%  find the correlation image for each landmark point:
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
NNeighbor = sum(adjM~=nLdmk+1,2);

% precompute the distance to center used in case of end points landmarks
hoodDist  =   abs(complex(X,Y));

% hoodAngle = angle(complex(X,Y));

% gaussian filter for the
hg = fspecial('gaussian', [3 3], 2);

%%
for iter= 1:20 %iterations
    ldmkMovingOld = ldmkMoving;
    for p = ldmkIdx
        if iter==2
            1;
        end
        hoodCoordTmp = bsxfun(@plus, hoodCoord, ldmkMoving(p,:));
        adjLdmkTmp = ldmkMoving(adjM(p,:),:);
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

%                     hoodCurv(:) = ...
%                         abs(sqrt(sum( bsxfun(@plus, -2*hoodCoordTmp, sum(adjLdmk(p,:,[1 2]),3)).^2 ,2)) - adjCurve(p,1)) ;

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

%                     hoodCurv(:) = ...
%                        (abs(sqrt(sum( bsxfun(@plus, -2*hoodCoordTmp, sum(adjLdmk(p,:,[1 2]),3)).^2 ,2)) - adjCurve(p,1))+ ...
%                         abs(sqrt(sum( bsxfun(@plus, -2*hoodCoordTmp, sum(adjLdmk(p,:,[1 3]),3)).^2 ,2)) - adjCurve(p,2))+ ...
%                         abs(sqrt(sum( bsxfun(@plus, -2*hoodCoordTmp, sum(adjLdmk(p,:,[2 3]),3)).^2 ,2)) - adjCurve(p,3))  )/3;

            otherwise
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
        
        %% plot 
        if true(1)%  && (p==12) % || p ==14)
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
                        ldmkMoving(1:end-1,:)  ).^2,2)))<1
        break
    end
end
    
    ldmkMoving = ldmkMoving(1:end-1,:);
    
    
end








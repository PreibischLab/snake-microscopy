


function [adjDist, adjCurve] = fAdjascentCalculations(ldmk, adjLdmk)

    iLdmk = true(size(ldmk,1),1);
    iLdmk(end) = false;

    % distance to adjascent points
    adjDist = sqrt( sum( bsxfun(@minus, adjLdmk, ldmk(iLdmk,:)).^2 ,2) );
    adjDist = permute(adjDist,[1 3 2]);

    adjCurve = nan(size(adjLdmk,1),3);
    % curvature
%     adjCurve(:,1) =  sqrt(sum( (sum(adjLdmk(:,:,[1 2]),3) - 2*ldmk(iLdmk,:)).^2 ,2));
%     adjCurve(:,2) =  sqrt(sum( (sum(adjLdmk(:,:,[1 3]),3) - 2*ldmk(iLdmk,:)).^2 ,2));
%     adjCurve(:,3) =  sqrt(sum( (sum(adjLdmk(:,:,[2 3]),3) - 2*ldmk(iLdmk,:)).^2 ,2));
%     adjCurve(isnan(adjCurve)) = 0;
    
    % oriented angle
    % (https://stackoverflow.com/questions/3486172/angle-between-3-points)
    V = bsxfun(@minus, adjLdmk(:,:,[1 2]), ldmk(iLdmk,:));
    vDot = sum(prod(V,3),2);
    vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
    adjCurve(:,1) = atan2(vCross, vDot);
    
    
    V = bsxfun(@minus, adjLdmk(:,:,[2 3]), ldmk(iLdmk,:));
    vDot = sum(prod(V,3),2);
    vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
    adjCurve(:,2) = atan2(vCross, vDot);
    
    adjCurve(:,3) = 0;
    
end

%%
% imagesc(xWindows, yWindows,template(yWindows, xWindows))
% 
% % imagesc(xWindows, yWindows,hoodCont)
% % imagesc(xWindows, yWindows,Y)
% axis equal tight
% 
% 
% %%
% hold on
% scatter(ldmkMoving(p,1), ldmkMoving(p,2),'+')
% line([ldmkMoving(p,1) adjLdmk(p,1,1)], [ldmkMoving(p,2) adjLdmk(p,2,1)])
% hold off
% 
% %%
% tmp = (-hoodSize):(hoodSize);
% 
% hooA = round(fRangeNorm(hoodAngle,0.51,100.49));
% u = unique(hooA);
% H = zeros(numel(tmp),numel(tmp),numel(u));
% for i=1:size(H,1)
%     for j=1:size(H,2)
%         H(i,j,:) = accumarray(hooA(:),reshape(template(yWindows, xWindows),[],1));
%     end
% end
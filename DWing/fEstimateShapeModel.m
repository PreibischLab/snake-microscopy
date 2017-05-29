% shapeModel = fEstimateShapeModel(landmark, adjIndices)
% shapeModel
% 
% 
% the function estimates :
% - the distance between the landmark and its adjascent points
%         (shapeModel.adjDistance)
% 
% - the oriented angle between  :
%           -  landmark and neighbor 1 and 2
%           -  landmark and neighbor 2 and 3
%       (shapeModel.adjCurvature)
% 
% and stores the indices of adjascent  points
%       (shapeModel.adjIndices)
% 
% paramters : 
%     - landmark: nLdmk x 2:
%             [x y coordinates of the landmarks
%     - adjMat:   nLdmk x 3:
%               index of the adjascent point for each landmark
%               or NaN on the last columns if less than 3 neighbors
% 

function shapeModel = fEstimateShapeModel(landmark, adjIndices)

nLdmk = size(landmark,1);

% create a false landmark [NaN, NaN] to make the calculations easy
landmark = [landmark; NaN NaN];

% refer to the false landmark for not existing adjascent points
adjIndices(isnan(adjIndices)) = nLdmk+1;


% table with dim1 = ldmk index,
%            dim2 = [X Y] coordinates 
%            dim3 = adj index,
% for easier calculations
adjLdmk = reshape(landmark(adjIndices,:),size(adjIndices,1),size(adjIndices,2),2);
adjLdmk = permute(adjLdmk,[1 3 2]);

% adjacency matrix
%     different adjascents
%    ¬
%   /-->   [x, y] coordinates
%   |
%   v
%   landmark index
% 
adjDistance = sqrt( sum( bsxfun(@minus, adjLdmk, landmark(1:end-1,:)).^2 ,2) );
adjDistance = permute(adjDistance,[1 3 2]);


% oriented angle between vectors between landmark and its adjascent points
% (https://stackoverflow.com/questions/3486172/angle-between-3-points)
adjCurvature = nan(size(adjLdmk,1),2);

% angle from [L->Adj1 and L->Adj2)
V = bsxfun(@minus, adjLdmk(:,:,[1 2]), landmark(1:end-1,:));
vDot = sum(prod(V,3),2);
vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
adjCurvature(:,1) = atan2(vCross, vDot);

% angle from [L->Adj2 and L->Adj3)
V = bsxfun(@minus, adjLdmk(:,:,[2 3]), landmark(1:end-1,:));
vDot = sum(prod(V,3),2);
vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
adjCurvature(:,2) = atan2(vCross, vDot);


shapeModel.adjDistance  = adjDistance;
shapeModel.adjCurvature = adjCurvature;
shapeModel.adjIndices   = adjIndices;


end
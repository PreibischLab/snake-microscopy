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
%       (shapeModel.adjCurveVal)
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

function shapeModel = fEstimateShapeModel(LandmarkIS, LdmkClass, WeightPerClass,...
                                          Segments, RelPosInSegments, ...
                                          AdjIdx, DistIdx, CurvIdx)
     
    nLdmk = size(LandmarkIS,1);


    shapeModel.segments         = Segments; 
    shapeModel.relPosInSegments = RelPosInSegments;

    % procruste superimposition
    modelIS = mean(LandmarkIS,3);
    if size(LandmarkIS,3) > 1
        modelIS = mean(LandmarkIS,3);
        [~,LandmarkBS] =  fProcrustesSupp(num2cell(LandmarkIS,[1 2]), false, modelIS);
        LandmarkBS = cell2mat(reshape(LandmarkBS,1,1,[]));
    else
        LandmarkBS = LandmarkIS;
    end
    
    shapeModel.startPosition = modelIS;
    

    % find the class of landmarks
    shapeModel.ldmkClass = LdmkClass;

    % weight per class
    crossPoints = shapeModel.ldmkClass==1;
    endPoints   = shapeModel.ldmkClass==2;
    semiPoints  = shapeModel.ldmkClass==3;



    % centroid distance
    modelBS = mean(LandmarkBS,3);
    shapeModel.centroidAvDist = ...
        mean(sqrt(sum(bsxfun(@minus, modelBS, mean(modelBS,1)).^2,2)));

    % adjascent angle (for semi landmarks)
    shapeModel.adjIndices = AdjIdx;
    [~, shapeModel.adjCurveVal, curvVar] = ...
        f_Curvature( LandmarkBS, AdjIdx, semiPoints );
    shapeModel.adjCurveWeight = f_VarianceToWeight(curvVar);
    
    % distance (for cross+end points)
    tmp = crossPoints | endPoints;
    if (~exist('DistIdx','var'))||isempty(DistIdx)
        % select points with less variance
        [DistIdx, distVar] = f_FindRobustNeighbors('distance', LandmarkBS, tmp, 2, shapeModel.centroidAvDist);
    end
    shapeModel.distIdx = DistIdx;
    [~, shapeModel.distSubVal, distSubVar] =...
        f_Distances(LandmarkBS, DistIdx, tmp, shapeModel.centroidAvDist);
  
    shapeModel.distSubWeight = f_VarianceToWeight(distSubVar);
        
    % curve (landmarks)
    if (~exist('CurvIdx','var'))||isempty(CurvIdx)
        % select points with less variance
        [CurvIdx, curvVar]= f_FindRobustNeighbors('curvature', LandmarkBS, tmp, 2, shapeModel.centroidAvDist); 
    end
    [~, shapeModel.curvSubVal, curvSubVar] = ...
            f_Curvature(LandmarkBS, CurvIdx, tmp);
    shapeModel.curvSubWeight = f_VarianceToWeight(curvSubVar);
    shapeModel.curvIdx = CurvIdx;

    
    shapeModel.weights = zeros(nLdmk,3);
    shapeModel.weights(crossPoints, :) = repmat(WeightPerClass(1,:), sum(crossPoints ), 1);
    shapeModel.weights(endPoints  , :) = repmat(WeightPerClass(2,:), sum(  endPoints ), 1);
    shapeModel.weights(semiPoints , :) = repmat(WeightPerClass(3,:), sum( semiPoints ), 1);
    if exist('distVar','var') && exist('curvVar','var')
        % change the weight for cross and end points
        i = reshape(find(tmp),1,[]);
        for idxL = i
        	t = f_VarianceToWeight({[sum(distVar{idxL}) sum(curvVar{idxL})]});
            t = t{1} .* shapeModel.weights(idxL,[1 2]);
            t = t /sum(t);
            s = sum(shapeModel.weights(idxL,[1 2]));
            shapeModel.weights(idxL,[1 2]) = t  *s;
        end
    end
    
    % process order
    if size(LandmarkIS,3)>1
        V = mean(sqrt(sum(bsxfun(@minus, LandmarkIS, modelIS).^2,2)),3);

        idx = find(ismember(LdmkClass, [1 2]));
        [~,i] = sort(V(idx));
        shapeModel.processOrder = idx(i);

        idx = find(LdmkClass==3);
        [~,i] = sort(V(idx));
        shapeModel.processOrder = [shapeModel.processOrder ; idx(i)];
    else
        shapeModel.processOrder = [ find(LdmkClass==1);...
                                    find(LdmkClass==2);...
                                    find(LdmkClass==3) ];
    end

end

function weight = f_VarianceToWeight(cVariance)

weight = cell(size(cVariance));
    for i=1:numel(weight)
        if numel(cVariance{i})>1
            weight{i} = 1-cVariance{i}/sum(cVariance{i});
        else
            weight{i} = 1;
        end
    end

end

function [idxSelected, Var] = f_FindRobustNeighbors(strOption, LandmarksBS, filter, keepN, centroidSize)
filter = reshape(filter,1,[]);
switch strOption
    case 'curvature'
        nLdmk = size(LandmarksBS,1);
        AdjIdx = cell(nLdmk,1);
        if islogical(filter)
            filter = reshape(find(filter),1,[]);
        end
        
        for l = filter
            AdjIdx{l} = nchoosek(setdiff(filter,l),2);
        end
        
        varianceWeight = f_VarianceWeightEstimate(mean(LandmarksBS,3), AdjIdx, filter);
        
        [~, ~, adjCurveVar] =  f_Curvature( LandmarksBS, AdjIdx, filter );

        idxSelected = cell(nLdmk,1);
        Var         = cell(nLdmk,1);
        for l = filter
            Var{l} = adjCurveVar{l}.*varianceWeight{l} * 3 / centroidSize;
            [~,i] = sort(Var{l});
            iSel = i(1:min(keepN, numel(i)));
            idxSelected{l} = AdjIdx{l}(iSel,:);
%             Var{l} = Var{l}(iSel);
        end
        
        
    case 'distance'
        nLdmk = size(LandmarksBS,1);
        DistIdx = cell(nLdmk,1);
        if islogical(filter)
            filter = find(filter);
        end
        
        for l = filter
            DistIdx{l} = setdiff(filter, l);
        end
        
        [~, ~, distSubVar] =  f_Distances(LandmarksBS, DistIdx, filter, centroidSize);
        
        idxSelected = cell(nLdmk,1);
        Var         = cell(nLdmk,1);
        for l = filter
            Var{l} = distSubVar{l};
            [~,i] = sort(distSubVar{l});
            iSel = i(1:min(keepN, numel(i)));
            idxSelected{l} = DistIdx{l}(iSel);
%             Var{l} = Var{l}(iSel);
        end
        
    otherwise
        error('%s not recognized as a valid option',strOption)
end

end


% function [curveValInd, curveValAv, curveWeight] = f_Curvature(landmarks, idx)
% 
% landmarks: array -> nLdmk x 2 x nIndividuals
%     landmarks coordinates
% idx: cell -> nLdmk x 1
%     index of neighbors
%     each cell element is nAngles x 2 (index of the landmarks)
% 
function [curveValInd, curveValAv, curveVariance] = f_Curvature(landmarks, idxAngle, idxSelect)
% idx cell: nLdmk x 1
%   each cell nAngle x 2 (2-> two points to which the central point is
%   compared)
    idxAngle = reshape(idxAngle,[],1);

    if islogical(idxSelect)
        idxSelect = find(idxSelect);
    end
    idxSelect = reshape(idxSelect,1,[]);
    
    %% calculate the different angles
    % do all calculation at once by constructing a big array of neighbors
    idxNeigh = cat(1, idxAngle{idxSelect});
    NNeighbor = cellfun(@(x) size(x,1), idxAngle(idxSelect));
    idxLandmark = arrayfun(@(x) x*ones(size(idxAngle{x},1),1),(1:size(idxAngle,1))','Unif',false);
    idxLandmark = cat(1,idxLandmark{idxSelect});
%     idxNeighbor = arrayfun(@(x) (1:size(idxAngle{x},1))',(1:size(idxAngle,1))','Unif',false);
%     idxNeighbor = cat(1,idxNeighbor{idxSelect});
    if size(idxNeigh,2) ~= 2
        error('angle is calculated between two neighbors')
    end
    curveValInd   = cell(size(landmarks,1), size(landmarks,3));
    
    for iInd=1:size(landmarks,3) % for each individual
        % table with dim1 = ldmk index,
        %            dim2 = [X Y] coordinates 
        %            dim3 = adj index,
        % for easier calculations
        % two neighbor and 2 coordinates (x and y)
        adjLdmk = reshape(landmarks(idxNeigh,:,iInd), sum(NNeighbor), 2,  2);
        adjLdmk = permute(adjLdmk,[1 3 2]); % put coordinate in second position
        
        % angle from [L->Adj1 and L->Adj2)
        V = bsxfun(@minus, adjLdmk(:,:,[1 2]), landmarks(idxLandmark,:,iInd));
        A = fAngleBtVectors(V);

        % reshape to get the values for each landmarks
        curveValInd(idxSelect,iInd) = mat2cell(A,NNeighbor);
%         distValInd(idxSelect,iInd)  = mat2cell(D,NNeighbor);        

    end
    curveValInd = cellfun(@(x) reshape(x,1,[]),curveValInd,'unif',false);
%     distValInd  = cellfun(@(x) reshape(x,1,[]),distValInd ,'unif',false);

    %% calculate average and variance for the different angles
    if nargout>1
        curveVariance = cell(size(landmarks,1),1);
        curveValAv  = cell(size(landmarks,1),1);
        if size(landmarks,3)>1
            n = size(landmarks,3)-1;

            for iL = idxSelect % for each landmark
                
                V = cat(1,curveValInd{iL,:});
                
                % angle difference with first individual
                A = fAngularDiff(V,V(1,:),'signed');
                A = bsxfun(@minus, A, mean(A,1));
                curveVariance{iL} = sqrt(sum(A.^2,1))/n;
                
                % convert angle to pixel units so that it is comparible to
                % the distance
                curveVariance{iL} = curveVariance{iL};
                
                % average angle: use the angle distance distance
                m = mean(A,1);
                curveValAv{iL} = V(1,:) + m;
            end
            
        else
            for iL = idxSelect % for each landmark
                curveValAv{iL}  = curveValInd{iL,1};
                curveVariance{iL} = ones(1,numel(curveValInd{iL,1}));
            end
            
        end
    end

end


function [distValInd, distValAv, distVariance] = f_Distances(landmarks, idxDist, idxSelect, centroidSize)
% idx cell nLdmk x 1
%   idx{iLdmk,iInd} nNeighbor{iLdmk,iInd} x 1

    idxDist = cellfun(@(x) reshape(x,[],1),idxDist,'unif',false);
    
    idxNeigh = cat(1,idxDist{idxSelect});
    idxNeighbor = cellfun(@(x) size(x,1),idxDist(idxSelect));
    idxLandmark = arrayfun(@(x) x*ones(size(idxDist{x},1),1),(1:size(idxDist,1))','Unif',false);
    idxLandmark = cat(1,idxLandmark{idxSelect});
    
    distValInd    = cell(size(landmarks,1), size(landmarks,3));
    for iInd=1:size(landmarks,3) % for every individual
        % table with dim1 = ldmk index,
        %            dim2 = [X Y] coordinates 
        %            dim3 = adj index,
        % for easier calculations
        adjLdmk = reshape(landmarks(idxNeigh,:,iInd),sum(idxNeighbor), 2);        
        
        D = sqrt( sum( (adjLdmk-landmarks(idxLandmark,:,iInd)).^2 ,2) );
        D = permute(D,[1 3 2])/centroidSize;
        
        % reshape to get the values for each landmarks
        distValInd(idxSelect,iInd) = mat2cell(D,idxNeighbor);
    end
    distValInd = cellfun(@(x) reshape(x,1,[]),distValInd,'unif',false);
    
    if nargout>1
        distValAv = cell(size(landmarks,1),1);
        distVariance = cell(size(landmarks,1),1);
        if size(landmarks, 3)>1 % if several individuals
            for iL=1:size(distValInd,1)% for each landmark
                if numel(idxDist{iL})>0
                    D = cat(1,distValInd{iL,:});
                    distValAv{iL}  = mean(D,1);
                    distVariance{iL} = std(D,0,1);
                end
            end
        else
            for iL=1:size(distValInd,1)% for each landmark
                distValAv{iL} = distValInd{iL,1};
                distVariance{iL} = ones(1,numel(distValInd{iL,1}));
            end
        end
    end


end

function varianceWeight = f_VarianceWeightEstimate(landmarks, idxAngle, idxSelect)
% estimate how much a variation in the position of the landmarks influence
% the curvature. This is used to let distance and angle being comparible

R = 10;
XYnei = R*[ cos((0:8)/4*pi)',...
            sin((0:8)/4*pi)' ];
XYnei = [0 0; XYnei];
varianceWeight = cell(size(landmarks,1),1);

    for iLdmk = idxSelect
        adjLdmk = reshape(landmarks(idxAngle{iLdmk},:,1), [], 2,  2);
        adjLdmk = permute(adjLdmk,[1 3 2]); % put coordinate in second position
        ldmkNei = bsxfun(@plus, XYnei, landmarks(iLdmk,:));
        
        neiIdx = arrayfun(@(x) x*ones(size(adjLdmk,1),1),1:size(ldmkNei,1),'Unif',false);
        neiIdx = cat(1,neiIdx{:});
        adjIdx = repmat((1:size(idxAngle{iLdmk},1))',size(XYnei,1),1);
        V = bsxfun(@minus, repmat(adjLdmk,10,1,1), ldmkNei(neiIdx,:));
        A = fAngleBtVectors(V);
        Atmp = zeros(1,size(idxAngle{iLdmk},1));
        %%
        for iA = 1:size(idxAngle{iLdmk},1)
            i = find(adjIdx==iA);
            Atmp(iA) = max(fAngularDiff(A(i),A(i(1)),'abs'));
        end
        %%
        varianceWeight{iLdmk} = R./Atmp;
    end

end



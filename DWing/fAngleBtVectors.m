function anglDiff = fAngleBtVectors(V)
% calulate the oriented angle bt these
% dot and cross product
vDot = sum(prod(V,3),2);
vCross = diff(V(:,[2 1],1).*V(:,:,2),1,2);
    
% get the angle
anglDiff = atan2(vCross, vDot);
end
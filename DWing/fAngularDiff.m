

function anglDiff = fAngularDiff(A1, A2, opt)
if (~exist('opt','var')) ||isempty(opt)
    opt = 'abs';
end
switch opt
    case 'abs'
        % warp it to that the results are comparible
        anglDiff = mod(bsxfun(@minus,A1 , A2), 2*pi);
        anglDiff = min(2*pi-anglDiff ,anglDiff );
        
    case 'signed'
        anglDiff = mod(bsxfun(@minus,A1 , A2)-pi,2*pi)-pi;
end

end
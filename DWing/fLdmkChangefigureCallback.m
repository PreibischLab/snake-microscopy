function fLdmkChangefigureCallback(hObject,callbackdata,...
                                   hSubs, iLdmk, iImg, hoodSize,...
                                   Istack, LdmkStack,...
                                   template, ldmk, force)
if ~exist('force','var')
    force = false;
end
boolRefresh = true(1);
persistent hPoint;
    iLdmk = evalin('base', 'iLdmk');
    iImg = evalin('base', 'iImg');
    LdmkStack = evalin('base', 'LdmkStack');

% below:
if ~force
    switch callbackdata.Key
        case 'f5'
            
        case 'leftarrow' 
            iImg = iImg-1;
            if iImg<1
                iImg = 1;
            end
        case 'rightarrow'
            iImg = iImg+1;
            if iImg>numel(Istack)
                iImg = numel(Istack);
            end
        case 'uparrow'
            iLdmk = iLdmk+1;
            if iLdmk>size(ldmk,1)
                iLdmk = size(ldmk,1);
            end
        case 'downarrow'
            iLdmk = iLdmk-1;
            if iLdmk<1
                iLdmk = 1;
            end
        case {'s' 'control' 'numpad0' 'insert'}
            if exist('hPoint','var') && (~isempty(hPoint))
                P = round(hPoint.getPosition());
                LdmkStack{iImg,2}(iLdmk,:) = round(P); 
                assignin('base', 'LdmkStack', LdmkStack);
                boolRefresh = true;
            end
        otherwise
            boolRefresh = false;
    end
end

if boolRefresh
    if ~isempty(hPoint)
        delete(hPoint)
    end
    rIm = (round(LdmkStack{iImg, 2}(iLdmk,2)-hoodSize)):...
          (round(LdmkStack{iImg, 2}(iLdmk,2)+hoodSize));
    cIm = (round(LdmkStack{iImg, 2}(iLdmk,1)-hoodSize)):...
          (round(LdmkStack{iImg, 2}(iLdmk,1)+hoodSize));

    rIm = rIm((rIm<=size(Istack{iImg},1)) & (rIm>0));
    cIm = cIm((cIm<=size(Istack{iImg},2)) & (cIm>0));

    rTp = (ldmk(iLdmk,2)-hoodSize):...
          (ldmk(iLdmk,2)+hoodSize);
    cTp = (ldmk(iLdmk,1)-hoodSize):...
          (ldmk(iLdmk,1)+hoodSize);

    rTp = rTp((rTp<=size(template,1)) & (rTp>0));
    cTp = cTp((cTp<=size(template,2)) & (cTp>0));

    axes(hSubs(1));
    imagesc(cTp,rTp, template(rTp,cTp));
    set(hSubs(1),'xtick',[])
    set(hSubs(1),'ytick',[])
    hold on
    scatter(ldmk(iLdmk,1),ldmk(iLdmk,2),'+','r');
    hold off
    axis equal tight

    axes(hSubs(2));
    imagesc(cIm,rIm, Istack{iImg}(rIm,cIm));
    set(hSubs(2),'xtick',[])
    set(hSubs(2),'ytick',[])
    hPoint = impoint(hSubs(2), LdmkStack{iImg, 2}(iLdmk,:));
    axis equal tight
    
    assignin('base', 'iLdmk', iLdmk);
    assignin('base', 'iImg', iImg);
    drawnow;
end

end
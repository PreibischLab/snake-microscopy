function fLdmkChangefigureCallback(hObject,callbackdata,...
                                   hSubs, iLdmk, iImg, hoodSize,...
                                   Istack, LdmkStack,...
                                   template, ldmk, force, ldmkLocalBestCell)
if (~exist('force','var')) || isempty(force)
    force = false;
end
boolRefresh = true(1);
persistent hPoint;
persistent plotBest;

iLdmk = evalin('base', 'iLdmk');
iImg = evalin('base', 'iImg');
LdmkStack = evalin('base', 'LdmkStack');
bDisplayBest = evalin('base', 'bDisplayBest');
if isempty(plotBest)
    plotBest = false;
end
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
            
        case 'space'
            plotBest = ~plotBest;
            
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
    scatter(ldmk(iLdmk,1),ldmk(iLdmk,2),'+','b');
    hold off
    axis equal tight

    axes(hSubs(2));
    imagesc(cIm,rIm, Istack{iImg}(rIm,cIm));
    set(hSubs(2),'xtick',[])
    set(hSubs(2),'ytick',[])
    axis equal tight
    if plotBest && (exist('ldmkLocalBestCell','var') || isempty(ldmkLocalBestCell))
        % shape
        h = zeros(size(bDisplayBest));
        strLegend = cell(size(h));
        if bDisplayBest(1)
            hold on
            [~,idx] = sort(ldmkLocalBestCell{iImg}(iLdmk).shape{2});
            idx = idx(1:min(2, numel(idx)));
            h(1) = scatter(ldmkLocalBestCell{iImg}(iLdmk).shape{1}(idx,1),...
                           ldmkLocalBestCell{iImg}(iLdmk).shape{1}(idx,2),300,'+','g');
            strLegend{1} = 'shape best position';
            hold off
        end
        
        % metric
        if bDisplayBest(2)
            map = [1 0 0; 1 0.5237 0; 1 0.8 0];
            hold on
            [~,idx] = sort(ldmkLocalBestCell{iImg}(iLdmk).metric{2});
            idx = idx(1:min(3, numel(idx)));
%             idx = MV<=median(MV);
            h(2) = scatter(ldmkLocalBestCell{iImg}(iLdmk).metric{1}(idx,1),...
                           ldmkLocalBestCell{iImg}(iLdmk).metric{1}(idx,2),300,map,...
                        '+');
            strLegend{2} = 'metric best position';
            hold off
        end
        
        if any(bDisplayBest)
            legend(h(bDisplayBest),strLegend)
        end
        
    end
    hPoint = impoint(hSubs(2), LdmkStack{iImg, 2}(iLdmk,:));
    hPoint.setColor([82 135 255]/255);

    
    assignin('base', 'iLdmk', iLdmk);
    assignin('base', 'iImg', iImg);
    drawnow;
end

end
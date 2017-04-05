

% path = ('/home/ella/Desktop/dataset');
path = ('/Users/ebahry/Desktop/dataset');

stacks = [2];

load('slicesNucs.mat');
sizSlicesNucs = size(slicesNucs,1);
for i=1:size(stacks,2)
    
	for j=1:80
        
        load('slicesNucsAll.mat');

%         img = imread([path '/EM_' num2str(stacks(i)) '_' num2str(stacks(i)+10) '.tiff'],j);
        img = imread([path '/full_worm_size10_part' num2str(stacks) '.tif'],j);
        startLoc = 4000;

        f=figure;imshow(img(:,startLoc:end)); hold on; impixelinfo;
        xlim([250 2000]);
        ylim([500 1200]);
        set(gcf, 'Position', [0 700 1800 1500]);
        
        slice = (stacks(i)-1)*80+j;
        title({'slice ', num2str(slice)});
        
        nucs = slicesNucs(:,slice);

        for nuc=1:size(nucs,1)
            temp = nucs{nuc};
            if ~isempty(temp)
                Xs = temp(:,1);
                Ys = temp(:,2);
                handl = plot(Ys-startLoc, Xs);
                text(Ys(1)-startLoc,Xs(1), num2str(nuc), 'Color','red','FontSize', 30);
            end
        end
        
        toEdit = input('edit?');
        while ~isempty(toEdit)
            
            [ytemp, xtemp] = ginput(4);
            nucNumOnNuc = isNucOnNuc([round(xtemp), round(ytemp) + startLoc], slice);
            
            Xs = nan;
            Ys = nan;
            h = imfreehand(); 
            pos = getPosition(h);
            Ys = round(pos(:,1) + startLoc);
            Xs = round(pos(:,2));
            slicesNucs{nucNumOnNuc,slice} = [Xs,Ys];

            kids = get(gca, 'Children');
            delete(kids(1:end-1));

            nucs = slicesNucs(:,slice);

            for nuc=1:size(nucs,1)
                temp = nucs{nuc};
                if ~isempty(temp)
                    Xs = temp(:,1);
                    Ys = temp(:,2);
                    handl = plot(Ys-startLoc, Xs);
                    text(Ys(1)-startLoc,Xs(1), num2str(nuc), 'Color','red','FontSize', 30);
                end
            end
            
            toEdit = input('edit?');
        end
            
        save('slicesNucsAll','slicesNucs');
        
        
%         delNucs = input('which to delete?');
%         for l=1:size(delNucs,2)
%             slicesNucs{delNucs(l),slice-isSN2} = [];
%         end
%         % Just rearraigning to have no empty cells in middle
%         a = slicesNucs(:,slice-isSN2); 
%         b=a(~cellfun('isempty',a));
%         sizB = size(b,1);
%         slicesNucs(1:sizB,slice-isSN2)=b;
%         c = cell(sizSlicesNucs-sizB,1);
%         slicesNucs(sizB+1:end,slice-isSN2) = c;
%         
%         save('slicesNucs2','slicesNucs');
%         
%         kids = get(gca, 'Children');
%         delete(kids(1:end-1));
%         
%         nucs = slicesNucs(:,slice-isSN2);
% 
%         for nuc=1:size(nucs,1)
%             temp = nucs{nuc};
%             if ~isempty(temp)
%                 Xs = temp(:,1);
%                 Ys = temp(:,2);
%                 handl = plot(Ys-startLoc, Xs);
%                 text(Ys(1)-startLoc,Xs(1), num2str(nuc), 'Color','red','FontSize', 30);
%             end
%         end
%         
%         addN = 1;
%         % How many new did we add
%         % newNucs holds all the new nucs we draw
%         newNucs = cell(50,1);
%         counterNew = 1;
%         
%         while ~isempty(addN)
%             addN = input('how many to add?');
% 
%             for l=1:addN
%                 Xs = nan;
%                 Ys = nan;
%                 h = imfreehand(); 
%                 pos = getPosition(h);
%                 Ys = round(pos(:,1) + startLoc);
%                 Xs = round(pos(:,2));
%                 
%                 newNucs{counterNew} = [Xs, Ys];
%                 counterNew = counterNew + 1;
% 
%                 slicesNucs{end-l+1,slice-isSN2} = [Xs,Ys];
%                 
%                 if (size(slicesNucs{end-l+1,slice-isSN2},1) < 10)
%                     slicesNucs{end-l+1,slice-isSN2} = [];
%                 else
%                     nucNumOnNuc = isNucOnNuc([Xs, Ys], slice);
% 
%                     if (~isnan(nucNumOnNuc))
%                     	slicesNucs{nucNumOnNuc,slice-isSN2} = [];
%                     end
%                 end
%             end
%             
%             % Just rearraigning to have no empty cells in middle
%             a = slicesNucs(:,slice-isSN2); 
%             b=a(~cellfun('isempty',a));
%             sizB = size(b,1);
%             slicesNucs(1:sizB,slice-isSN2)=b;
%             c = cell(sizSlicesNucs-sizB,1);
%             slicesNucs(sizB+1:end,slice-isSN2) = c;
% 
%             save('slicesNucs2','slicesNucs');
%             
%             kids = get(gca, 'Children');
%             delete(kids(1:end-1));
% 
%             nucs = slicesNucs(:,slice-isSN2);
% 
%             for nuc=1:size(nucs,1)
%                 temp = nucs{nuc};
%                 if ~isempty(temp)
%                     Xs = temp(:,1);
%                     Ys = temp(:,2);
%                     handl = plot(Ys-startLoc, Xs);
%                     text(Ys(1)-startLoc,Xs(1), num2str(nuc), 'Color','red','FontSize', 30);
%                 end
%             end
%         end
%         
%         if (j<0)
%             c = newNucs(~cellfun('isempty',newNucs));
%             for k = 1:size(c,1)
%                 % See if we have the new drawn nuc in the next annotated slice 
%                 nucNumInAnnot = isNucOnNuc(newNucs{k}, stacks(i)+5);
%                 % If so, see if we have it in the between slices and if not interpolate
%                 if ~isnan(nucNumInAnnot)
%                     % then iterpolate:
%                     % This function gets the slice number we're on, the
%                     % next annotation slice number, and the 2 nucs (on those 2 slices (so its the same nuc in 2 slices)) we interpolate between
%                     interpolNucBetSlices(slice, stacks(i)+5, newNucs{k}, slicesNucs{nucNumInAnnot,stacks(i)+5-isSN2})
%                 end
%             end
%         end
        close all;
	end
   
end


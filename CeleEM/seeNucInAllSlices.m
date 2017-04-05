close all;
clear;
load('OldSlicesNucsAndNucsNum')
uniqNucsNum = unique(nucsNum);
uniqNucsNum = uniqNucsNum(uniqNucsNum>0);

path = ('/Users/ebahry/Desktop/dataset');
%path = ('/home/ella/Desktop/dataset');

count = 0;
for i=1:size(uniqNucsNum)
    % i is nut my nuc number num is nucNum
    num = uniqNucsNum(i)
    
    [row, column] = find(nucsNum == num); 
    
    for j=1:size(row,1)
        
        nuc = slicesNucs{row(j), column(j)};
        
        part = ceil((column(j)+0)/80);
        slice = mod(column(j)+0,80);
        if (slice==0)
           slice = 80;
        end
        
        if j==1
            showFirstNuc = [nuc(1,1) nuc(1,2)];
        end
        
        addBuffer = 200;
        img = imread([path '/full_worm_size10_part' num2str(part) '.tif'], slice);
        f=figure;imshow(img(nuc(1,1)-addBuffer:nuc(1,1)+addBuffer,nuc(1,2)-addBuffer:nuc(1,2)+addBuffer)); hold on;
        %plot(nuc(1,2)-nuc(1,2)+addBuffer, nuc(1,1)-nuc(1,1)+addBuffer, '*');
        plot(showFirstNuc(2)-nuc(1,2)+addBuffer, showFirstNuc(1)-nuc(1,1)+addBuffer, '*');
        plot(nuc(:,2)-nuc(1,2)+addBuffer, nuc(:,1)-nuc(1,1)+addBuffer, '*');
        set(gcf, 'Position', [0 900 600 600]);
% %         pause();

    end
    
    pause();
    
    close all;
end


% %% to delete a nuc by its number:
% clear;
% load('nucsNumNew');
% load('slicesNucs');
% b = nucsNum==% NUC YOU WANT TO DELETE;
% nucsNum(b==1) = 0;
% c = sum(b(:));
% slicesNucs(b==1) = cell(1,c);
% save('slicesNucs','slicesNucs');
% save('nucsNum','nucsNum');

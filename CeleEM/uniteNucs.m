% % This script goes through the location of all nucs in each slice and
% % and unites the nuc from different slices 
% % by giving each nuc a unique number in "nucsNum"
% % After that it creates nucs4render where each nuc is a row
% % with all z x y of each nuc 
% 
% close all;
% clear;
%
% 
% % Init a matrix in which every nuc in slicesNucs will get a number
% % So we can trace the nuc through the slices
% nucsNum = zeros(maxNnucsInSlice,slicesN);
% 
% currNuc = 1;
% 
% for i=1:slicesN
%     tic
%     for j=1:maxNnucsInSlice
%         if (nucsNum(j,i)==0)
%            if ~(isempty(slicesNucs{j,i}))
%                
%                nucsNum(j,i) = currNuc;
%                
%                flag = 1;
%                slicesMore = 1;
%                
%                currContour = slicesNucs{j,i};
% 
%                while (flag~=0)
% 
%                    nextSlice = i+slicesMore;
%                    
%                    nucNumInOtherSlice = isNucOnNuc(currContour, nextSlice);
%                    
%                    if (isnan(nucNumInOtherSlice))
%                       flag = 0;
%                       currNuc = currNuc+1;
%                    else
%                        nucsNum(nucNumInOtherSlice(1), nextSlice) = currNuc;
%                        currContour = slicesNucs{nucNumInOtherSlice(1), nextSlice};
%                        slicesMore = slicesMore+1;
%                    end
%                end  
%            end
%         end
%     end
%     
%     toc
%     if (mod(i,10)==0)
%         save(['nucsNum' num2str(i)],'nucsNum');
%     end
%     
% end
% save('nucsNum','nucsNum');

load('slicesNucs');
load('nucsNum');

sizNucsNum = size(nucsNum);

nucsN = max(nucsNum(:));

nucs4render = cell(nucsN,30000);

for i=1:sizNucsNum(1)
    for j=1:sizNucsNum(2)
        
        if (nucsNum(i,j)~=0)
           row = nucsNum(i,j);
           
           nucInSlice = slicesNucs{i,j};
           howMany2insert = size(nucInSlice,1);
           where2insert = find(cellfun(@isempty,nucs4render(row,:)),1);
               
           for k=1:howMany2insert
               nucs4render{row,where2insert+k-1} = [j,nucInSlice(k,1),nucInSlice(k,2)];
           end
            
        end
        
    end
end

save('nucs4render','nucs4render')
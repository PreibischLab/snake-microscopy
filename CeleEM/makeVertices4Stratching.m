load('xyzStart.mat');
load('nucsContours4render');

fileID = fopen(['nucs4alignment.txt'],'w');

for i=1:721
    
    nucsContours = nucsContours4render(i,:);
    numOfSices = size(find(~cellfun(@isempty,nucsContours)),2);

    for j=1:numOfSices
        
        sizJ = size(nucsContours{j},2);
        
        for k=1:sizJ
            v = nucsContours{j}{k};
            siz = size(v,1);

            for l=1:siz
                fwrite(fileID,['Nuc' num2str(i) ' ' num2str(v(l,1) + xyzStart(i,1)) ' ' num2str(v(l,2) + xyzStart(i,2)) ' ' num2str((j+xyzStart(i,3))*2.25) char(10)]);
            end
        end

    end
end

fclose(fileID);
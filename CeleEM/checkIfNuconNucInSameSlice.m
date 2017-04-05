clear;
load('nucsNum.mat');
load('slicesNucs.mat');

for i=501:580
    i
    sliceNucs = slicesNucs(:,i);
    nucsN = sum(~cellfun(@isempty,sliceNucs));
    
    for j=1:nucsN
       
        nuc = sliceNucs{j};
        nucNumInOtherSlice = isNucOnNuc(nuc, i);
        
        if nucNumInOtherSlice ~= j
            
            if size(unique(sliceNucs{j},'rows'),1)<11
                slicesNucs{j,i}=[];
                nucsNum(j,i)=0;
                save('slicesNucs','slicesNucs');
                save('nucsNum','nucsNum');
                deleteEmpty;
                load('nucsNum.mat');
                load('slicesNucs.mat');
                display('fuuuuck!!!!!!!!!');
            else
                nucNumInOtherSlice
                load('nucsNum.mat');
                load('slicesNucs.mat');
            end
        end    
    end
end
    

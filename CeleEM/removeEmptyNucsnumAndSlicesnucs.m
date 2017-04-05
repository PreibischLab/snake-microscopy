load('OldSlicesNucsAndNucsNum');

for i=500:580

    a = slicesNucs(:,i); 
    b=a(~cellfun('isempty',a));
    sizB = size(b,1);
    slicesNucs(1:sizB,i)=b;
    c = cell(150-sizB,1);
    slicesNucs(sizB+1:end,i) = c;
    
    f = nucsNum(:,i);
    [loc, loc1, g] = find(f);
    h = f(f==0);
    p = [g;h];
    nucsNum(:,i) = p;
            
end

save('OldSlicesNucsAndNucsNum','slicesNucs','nucsNum');
        
        
files = dir('finalNuc*');
files = {files().name};
for i=1:184
    
    load(files{i});
    nuc(:,1) = nuc(:,1)/100;
    nuc(:,2) = nuc(:,2)/100;
    nuc(:,3) = nuc(:,3)/10;
    save(['a' files{i}],'nuc');
    
    name = files{i};
    name = name(1:end-4);
    csvwrite(['aa' name '.csv'], nuc);
    
end
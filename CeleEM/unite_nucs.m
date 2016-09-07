clear;

for i=101:184
    
    files = dir(['nuc' num2str(i) '_*']);
    files = {files().name};
    
    load(files{1});
    name = files{1};
    startSub = strfind(name,'e') + 1;
    endSub = strfind(name,'.') - 1;
    subStr = str2num(name(startSub:endSub));
    
    siz = size(Xs,1);
    
    Xss = [Xs(2:end); Xs(2)];
    Yss = [Ys(2:end); Ys(2)];
    Zss = nan;
    Zss(1:siz,1) = subStr;
    
    
    for j= 2:size(files,2)-1;
        load(files{j});
        name = files{j};
        startSub = strfind(name,'e') + 1;
        endSub = strfind(name,'.') - 1;
        subStr = str2num(name(startSub:endSub));

        siz = size(Xs,1);

        Xss = [Xss;Xs(2:end);Xs(2)];
        Yss = [Yss;Ys(2:end);Ys(2)];
        Zs = nan;
        Zs(1:siz,1) = subStr;
        Zss = [Zss;Zs];
              
    end
    
    startSub = strfind(name,'c') + 1;
    endSub = strfind(name,'_') - 1;
    subStr = str2num(name(startSub:endSub));
    
    nuc = [Xss,Yss,Zss];
    save(['finalNuc' num2str(subStr)], 'nuc');
    
end

    
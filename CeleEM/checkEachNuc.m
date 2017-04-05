% % Inspects each nucleous through all it's slices and checks if needs editing

path = ('/Users/ebahry/Desktop/dataset');
path = ('/home/ella/Desktop/dataset');

load('nucsNum');
load('slicesNucs');

maxNuc = max(nucsNum(:));

for i=1:max(nucsNum(:))
   
    [row, column] = find(nucsNum == i); 
    
    if (~isempty(row))
    
        [firstSlice, indFS] = min(column);
        [lastSlice, indLS] = max(column);

        nuc = slicesNucs{row(indLS), lastSlice};

        part = ceil(lastSlice/80);
        slice = mod(lastSlice,80);
        if (slice==0)
            slice = 80;
        end

        [i, part, slice]

        img = imread([path '/full_worm_size10_part' num2str(part) '.tif'], slice);
        f=figure;imshow(img(nuc(1,1)-200:nuc(1,1)+200,nuc(1,2)-200:nuc(1,2)+200)); hold on;
        plot(nuc(1,2)-nuc(1,2)+200, nuc(1,1)-nuc(1,1)+200, '*');
        set(gcf, 'Position', [0 900 600 600]);

        pause(0.01);
        showNext = 1;
        while (showNext)
            
            nextOnNuc = 0;

            if (slice==80)
                part = part+1;
                slice = 1;
            else
                slice = slice + 1;
            end
            
            windowSize = 200;

            img = imread([path '/full_worm_size10_part' num2str(part) '.tif'], slice);
            f=figure;imshow(img(nuc(1,1)-200:nuc(1,1)+200,nuc(1,2)-windowSize:nuc(1,2)+windowSize)); hold on;
            set(gcf, 'Position', [900 900 600 600]);   

            toEdit = input('to edit?');

            if (toEdit)
                Xs = nan;
                Ys = nan;
                h = imfreehand(); 
                pos = getPosition(h);
                Ys = round(pos(:,1)) + nuc(1,2)-windowSize;
                Xs = round(pos(:,2)) + nuc(1,1)-200;

                sliceInFull = (part-1)*80+slice;
                isThereANuc = isNucOnNuc([Xs,Ys], sliceInFull);

                if (isnan(isThereANuc))

                    emptyCell = find(cellfun(@isempty,slicesNucs(:,sliceInFull)),1);
                    slicesNucs{emptyCell, sliceInFull} = [Xs,Ys];
                    nucsNum(emptyCell, sliceInFull) = i;

                else
                    nucInSliceNum = nucsNum(isThereANuc(1), sliceInFull);
                    
                    nucsNum(nucsNum==nucInSliceNum(1)) = i;
                    nextOnNuc = 1;

                    [row, column] = find(nucsNum == i); 
                    [firstSlice, indFS] = min(column);
                    [lastSlice, indLS] = max(column);
                    nuc = slicesNucs{row(indLS), lastSlice};
                    part = ceil(lastSlice/80);
                    slice = mod(lastSlice,80);

                end
                save('slicesNucsAll','slicesNucs');
                save('nucsNumNew','nucsNum');
            end
            
            if (~nextOnNuc)
                showNext = input('show next?');
                set(gcf, 'Position', [0 900 600 600]);
            else
                showNext = 1;
            end
            
        end

        close all;
    
    end
    
end
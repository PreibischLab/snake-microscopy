load('nucsContours4render');
load('xyzStart.mat');

for i=513:513
    
    Vs = [1,1,1];

    if nucsWithDoubleContour(i) == false
        nucsContours = nucsContours4render(i,:);
        numOfSices = size(find(~cellfun(@isempty,nucsContours)),2);

        % Running over all slices to extract all vertices:
        for j=1:numOfSices
            % contour of j slice:
            vSlice = nucsContours{j};
            % Adding the y axis:
            vSlice(:,3) = j;

            Vs = [Vs; vSlice];

        end
        Vs(1,:) = [];
        Vs(:,4) = 1:size(Vs,1);
        
        Fs = [1,1,1];
        
        thisSlice = Vs(Vs(:,3)==1,:);
        for j = 2:size(thisSlice,1)-1
            Fs(end+1,:) = [thisSlice(1,4), thisSlice(j,4), thisSlice(j+1,4)];
        end
        thisSlice = Vs(Vs(:,3)==numOfSices,:);
        for j = 2:size(thisSlice,1)-1
            Fs(end+1,:) = [thisSlice(1,4), thisSlice(j,4), thisSlice(j+1,4)];
        end
        
        for j=1:numOfSices-1
            thisSlice = Vs(Vs(:,3)==j,:);
            nextSlice = Vs(Vs(:,3)==j+1,:);
            
            % If the two contours have the same amount of points:
            if size(thisSlice,1)==size(nextSlice,1)
                for k=1:size(thisSlice,1)-1
                    Fs(end+1,:) = [thisSlice(k,4), nextSlice(k,4), nextSlice(k+1,4)];
                    Fs(end+1,:) = [thisSlice(k,4), thisSlice(k+1,4), nextSlice(k+1,4)];
                end
                Fs(end+1,:) = [thisSlice(k+1,4), nextSlice(k+1,4), nextSlice(1,4)];
                Fs(end+1,:) = [thisSlice(k+1,4), thisSlice(1,4), nextSlice(1,4)];    
            % If the current contour has twice as many points than next:
            elseif ((size(thisSlice,1)/2)==size(nextSlice,1) || (size(thisSlice,1))==size(nextSlice,1)/2)
                if (size(thisSlice,1)/2)==size(nextSlice,1)
                    temp = thisSlice;
                    thisSlice = nextSlice;
                    nextSlice = temp;
                end
                for k=1:size(thisSlice,1)-1
                    Fs(end+1,:) = [thisSlice(k,4), nextSlice(2*k-1,4), nextSlice(2*k,4)];
                    Fs(end+1,:) = [thisSlice(k,4), nextSlice(2*k,4), nextSlice(2*k+1,4)];
                    Fs(end+1,:) = [thisSlice(k,4), thisSlice(k+1,4), nextSlice(2*k+1,4)];
                end
                Fs(end+1,:) = [thisSlice(k,4), nextSlice(2*k-1,4), nextSlice(2*k,4)];
                Fs(end+1,:) = [thisSlice(k,4), nextSlice(2*k,4), nextSlice(1,4)];
                Fs(end+1,:) = [thisSlice(k,4), thisSlice(1,4), nextSlice(1,4)];
                
            elseif ((size(thisSlice,1)/4)==size(nextSlice,1) || (size(thisSlice,1)==size(nextSlice,1)/4))
                if (size(thisSlice,1)/4)==size(nextSlice,1)
                    temp = thisSlice;
                    thisSlice = nextSlice;
                    nextSlice = temp;
                end
                for k=1:size(thisSlice,1)-1
                    Fs(end+1,:) = [thisSlice(k,4), nextSlice(3*k-2,4), nextSlice(3*k-1,4)];
                    Fs(end+1,:) = [thisSlice(k,4), nextSlice(3*k-1,4), nextSlice(3*k,4)];
                    Fs(end+1,:) = [thisSlice(k,4), nextSlice(3*k,4), nextSlice(3*k+1,4)];
                    Fs(end+1,:) = [thisSlice(k,4), thisSlice(k+1,4), nextSlice(3*k+1,4)];
                end
                Fs(end+1,:) = [thisSlice(k,4), nextSlice(3*k-2,4), nextSlice(3*k-1,4)];
                Fs(end+1,:) = [thisSlice(k,4), nextSlice(3*k-1,4), nextSlice(3*k,4)];
                Fs(end+1,:) = [thisSlice(k,4), nextSlice(3*k,4), nextSlice(1,4)];
                Fs(end+1,:) = [thisSlice(k,4), thisSlice(1,4), nextSlice(1,4)];
             elseif ((size(thisSlice,1)/8)==size(nextSlice,1) || (size(thisSlice,1)==size(nextSlice,1)/8))
                if (size(thisSlice,1)/4)==size(nextSlice,1)
                    temp = thisSlice;
                    thisSlice = nextSlice;
                    nextSlice = temp;
                end
                for k=1:size(thisSlice,1)-1
                    Fs(end+1,:) = [thisSlice(k,4), nextSlice(4*k-3,4), nextSlice(4*k-2,4)];
                    Fs(end+1,:) = [thisSlice(k,4), nextSlice(4*k-2,4), nextSlice(4*k-1,4)];
                    Fs(end+1,:) = [thisSlice(k,4), nextSlice(4*k-1,4), nextSlice(4*k,4)];
                    Fs(end+1,:) = [thisSlice(k,4), nextSlice(4*k,4), nextSlice(4*k+1,4)];
                    Fs(end+1,:) = [thisSlice(k,4), thisSlice(k+1,4), nextSlice(4*k+1,4)];
                end
                Fs(end+1,:) = [thisSlice(k,4), nextSlice(4*k-3,4), nextSlice(4*k-2,4)];
                Fs(end+1,:) = [thisSlice(k,4), nextSlice(4*k-2,4), nextSlice(4*k-1,4)];
                Fs(end+1,:) = [thisSlice(k,4), nextSlice(4*k-1,4), nextSlice(4*k,4)];
                Fs(end+1,:) = [thisSlice(k,4), nextSlice(4*k,4), nextSlice(1,4)];
                Fs(end+1,:) = [thisSlice(k,4), thisSlice(1,4), nextSlice(1,4)];
            else
                [i, j, size(thisSlice,1), size(nextSlice,1)]
            end
            
        end
        
        Fs(1,:) = [];
        
        fileID = fopen(['nuc' num2str(i) '.obj'],'w');
        for j=1:size(Vs,1)
            fwrite(fileID,['v ' num2str(Vs(j,1)+xyzStart(i,1)) ' ' num2str(Vs(j,2)+xyzStart(i,2)) ' ' num2str((Vs(j,3)+xyzStart(i,3))*2.25) ' 1.0' char(10)]);
        end

        for j=1:size(Fs,1)
            fwrite(fileID,['f ' num2str(Fs(j,1)) ' ' num2str(Fs(j,2)) ' ' num2str(Fs(j,3)) ' 1.0' char(10)]);
        end
        
        fclose(fileID);

    end
end



% This function checks the masks - 
% Checks if there is a slice that the nuc is in the wrong location:
% If so deletes it.
% For those nucs the imagej unite nuc should be ran again.
% Checks for holes.

path = '/Users/ebahry/Desktop/EM/nucsMasks/';
imgName = 'combinedNuc';
imgType = '.tif';

pathImg2Delete = '/Users/ebahry/Desktop/EM/EMnucsMasksJPG/';
emptyImg = imread([pathImg2Delete 'empty.jpg']);
arrayOfCorrected = zeros(50,1);
counterCorr = 1;
counter = 1;
for i=513:513
    fileName = ([path imgName num2str(i) imgType]);
    
    stackSize = numel(imfinfo(fileName));
    
    im1 = imread(fileName, 1);
    im1 = im2bw(im1, 0.5);
    
    for j=2:stackSize
        if counter>1
            counter = counter-1;
        else
            im2 = imread(fileName, j);
            im2 = im2bw(im2, 0.5);
            if max(max(im1+im2))~=2
                    
    %             % Check if this nuc has a slice that doesn't match the location
    %             % of the rest of the nuc:
    %             if (max(max(im1))==1 && max(max(im2==1)))
    %                 if j==2
    %                     i
    %                     j
    %                 else
    %                     im2 = emptyImg;
    %                     im2 = im2bw(im2, 0.5);
    %                     imwrite(im2, [pathImg2Delete '/nuc' num2str(i) '_' num2str(j) '.jpg']);
    %                     arrayOfCorrected(counter) = i;
    %                     counter = counter + 1;
    %                  end
    %             end
                counter = 1;
                bool = 1;
                dont = 0;
                while (bool==1)
                    if ~(stackSize>j+counter-1)
                        bool = 0;
                        dont = 1;
                    else
                        imTemp = imread(fileName, j+counter);
                        imTemp = im2bw(imTemp, 0.5);
                        if max(max(imTemp))==0
                            counter = counter+1;
                        else
                            bool = 0;
                        end 
                    end
                end
                if (j==2) && (max(max(im1))==0)
                else
                    if (dont == 0)
                        % find contour in first image
                        [tempX, tempY] = find(im1,1); 
                        edgeIm1 = bwtraceboundary(im1, [tempX tempY],'N');
                        edgeIm1 = edgeIm1(1:end-1,:);
                        % find contour in second image
                        [tempX, tempY] = find(imTemp,1); 
                        edgeImTemp = bwtraceboundary(imTemp, [tempX tempY],'N');
                        edgeImTemp = edgeImTemp(1:end-1,:);

                        s = RandStream('mt19937ar','Seed',0);
                        if size(edgeIm1,1)>size(edgeImTemp,1)
                            diffe = size(edgeIm1,1)-size(edgeImTemp,1);
                            r = randperm(s, size(edgeIm1,1), diffe);
                            edgeIm1(r,:)=[];
                        else
                            diffe = size(edgeImTemp,1)-size(edgeIm1,1);
                            r = randperm(s, size(edgeImTemp,1), diffe);
                            edgeImTemp(r,:)=[];
                        end

                        newXs = nan(1,counter);
                        newYs = nan(1,counter);
                        for k=1:size(edgeImTemp,1)
                            newXs(k,:) = linspace(edgeIm1(k,1), edgeImTemp(k,1), counter);
                            newYs(k,:) = linspace(edgeIm1(k,2), edgeImTemp(k,2), counter);
                        end
                        for k=1:counter
                            polyg = poly2mask(newYs(:,k), newXs(:,k), size(im1,1), size(im1,2));
                            imwrite(polyg, [pathImg2Delete '/nuc' num2str(i) '_' num2str(j+k-1) '.jpg']);
                        end
                        im2 = polyg;
                        arrayOfCorrected(counterCorr) = i;
                        counterCorr = counterCorr + 1;
                    end
                end
            end
            im1 = im2;
            
        end
    end
    
end
% arrayOfCorrected
save('arrayOfCorrected','arrayOfCorrected');

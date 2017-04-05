pathRead = '/Users/ebahry/Desktop/EM/nucsMasks/';
fileName = 'combinedNuc';
fileFormat = '.tif';

minPoints = 10;
maxContour4minPoints = 50;

nucsContours4render = cell(720, 100);
nucsWithDoubleContour = false(720,1);

for i=1:721

    imfile = ([pathRead fileName num2str(i) fileFormat]);
    imInfo = imfinfo(imfile);
    zSize = numel(imInfo);

    for j=1:zSize
        img = imread(imfile, j);
        img = imbinarize(img);

        [temp, countContours]=bwlabel(img);

        while (sum(sum(temp==1))<15 || (sum(sum(temp==2))<15 && (any(any(temp==2))))  || sum(sum(temp==3))<15 && (any(any(temp==3))))
            if sum(sum(temp==1))<15
                img(temp==1) = false;
            end
            if sum(sum(temp==2))<15
                img(temp==2) = false;
            end
            if sum(sum(temp==3))<15
                img(temp==3) = false;
            end
            [temp, countContours]=bwlabel(img);
        end

        count = 1;
        while (any(any(temp==count)))
            img1 = img;
            img1(temp~=count) = false;
            [tempX, tempY] = find(img1,1); 
            %figure; imshow(img); 
            edgePolyg = bwtraceboundary(img1, [tempX tempY],'N');
            % Delete the last point - because its equal to first point.
            edgePolyg = edgePolyg(1:end-1,:);

            sizContour = size(edgePolyg,1);

            if (sizContour<50)
                num = 0;
            elseif (sizContour<100)
                num = 1;
            elseif (sizContour<200)
                num = 2;
            elseif (sizContour<400)
                num = 3;
            elseif (sizContour<800)
                num = 4;
            end


            sample = (sizContour/(2^num*minPoints));
            % return the contour:
            Xs = edgePolyg(1:sample:end, 1);
            Ys = edgePolyg(1:sample:end, 2);

            %figure; imshow(img); hold on; plot(Ys(1),Xs(1),'o');
            %pause();

            nucsContours4render{i,j}{count} = [Xs, Ys];

            if (count>1)
                nucsWithDoubleContour(i) = true;
            end

            count = count + 1;

        end
    end

end
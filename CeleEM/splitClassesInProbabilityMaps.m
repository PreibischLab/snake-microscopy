% Splitting the probablity maps to the different classes 

function splitClassesInProbabilityMaps(imgName)

    for i=1:22
       img = imread(imgName,i);
       mat(:,:,i) = img;
    end
    mat = double(mat);
    imwrite(mat(:,:,1), '88probability_class1.tif');
    imwrite(mat(:,:,2), '88probability_class2.tif');
    for i=3:2:22
        imwrite(mat(:,:,i), '88probability_class1.tif', 'writemode', 'append');
    end
    for i=2:2:22
        imwrite(mat(:,:,i), '88probability_class2.tif', 'writemode', 'append');
    end
end
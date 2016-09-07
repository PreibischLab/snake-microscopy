% Splitting the probablity maps to the different classes 

function splitClassesInProbabilityMaps(imgName)

    for i=1:36
       img = imread(imgName,i);
       mat(:,:,i) = img;
    end
    mat = double(mat);
    imwrite(mat(:,:,1), 'probability_class1.tif');
    imwrite(mat(:,:,2), 'probability_class2.tif');
    for i=3:2:36
        imwrite(mat(:,:,i), 'probability_class1.tif', 'writemode', 'append');
    end
    for i=2:2:36
        imwrite(mat(:,:,i), 'probability_class2.tif', 'writemode', 'append');
    end
end
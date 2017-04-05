% this function first gets rid of redundent points and nucs that are too small:
% then resampls the contour every "sample" point.

load('slicesNucs');
load('nucsNum');
path = '/Users/ebahry/Desktop/dataset/';
path = ('/home/ella/Desktop/dataset/');
files = dir(path);
files = {files().name};
file = files{3};
img = imread([path file]);
imgSize = size(img);

sample = 2;

% First getting rid of redundent points and nucs that are too small:
for i=1:size(slicesNucs,1)
    for j=1:size(slicesNucs,2)
        
        if (~isempty(slicesNucs{i,j}))
            % remove redundent rows:
            temp = slicesNucs{i,j};
            temp = unique(temp,'rows','stable');
            slicesNucs{i,j} = temp;
            
            % Check if nuc is too small - then deletes it.
            if (size(slicesNucs{i,j},1) < 5)
                slicesNucs{i,j} = [];
                nucsNum(i,j) = 0;
            else
                % Resemples the contour:
                temp = resampleContour(imgSize, temp(:,1), temp(:,2), sample);
                temp = unique(temp,'rows','stable');
                slicesNucs{i,j} = temp;
                
            end
        end
    end
end

save('slicesNuc','slicesNucs');
save('nucsNum', 'nucsNum');
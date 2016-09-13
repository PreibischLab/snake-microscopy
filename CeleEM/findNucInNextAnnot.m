clc;    % Clear the command window.
close all;  % Close all figures
clear; % Clear workspace variables.

path = ('/home/ella/Desktop/dataset');

% % Convert all annotations to contours 
% contourSample = 4;
% for annotSlice=108:10:198
%     
%     labelImg = imread([path, '/1Label_size10_slice', num2str(annotSlice), '.tif']);
% 
%     nucsContours = label2contour(labelImg, contourSample);
%     save(['slice', num2str(annotSlice), '_Nucs'], 'nucsContours');
% end

% Set annotated slice number:
annotSlice = 28;

% Show image of next annotated slice
imgNext = imread([path, '/EM_', num2str(annotSlice+10), '.tif']);
gcf1 = figure(1); imshow(imgNext); hold on;
f1x = gca;

% Plot all annotated nucs on slice image
load(['slice' num2str(annotSlice+10) '_Nucs']);
for i=1:size(nucsContours,1)
    temp = nucsContours{i};
    Xs = temp(:,1);
    Ys = temp(:,2);
    plot(Ys, Xs);
    text(Ys(1),Xs(1), num2str(i), 'Color','red','FontSize', 30);
end

pause(0.5);

% Show image of the slice between the two annotated slices
imgMiddle = imread([path, '/EM_', num2str(annotSlice-5), '_', num2str(annotSlice+5) '.tiff'], 11);
gcf2 = figure(2); imshow(imgMiddle); hold on;
f2x = gca;

pause(0.5);

% Show image of annotated slice
imgPrev = imread([path, '/EM_', num2str(annotSlice), '.tif']);
gcf3 = figure(3); imshow(imgPrev); hold on;

load(['slice' num2str(annotSlice) '_Nucs']);
% Init cell array for mucs in mid slice
midNucs = cell(size(nucsContours,1),1);
% Init array for next slice nucs matching:
nextNucs = nan(size(nucsContours,1),1);
for i=1:size(nucsContours,1)
    % Plot all annotated nucs on slice image
    temp = nucsContours{i};
    Xs = temp(:,1);
    Ys = temp(:,2);
    plot(Ys, Xs);
    text(Ys(1),Xs(1), num2str(i), 'Color','red','FontSize', 30);

     % Show current nuc on annotated slice
    loc = round(size(Xs,1)/4);
    xlim([Ys(loc)-100 Ys(loc)+100]);
    ylim([Xs(loc)-100 Xs(loc)+100]);
    set(gcf3, 'Position', [0 0 600 600]); % Maximize figure
    % Show same location on next annotated slice
    f1x.XLim = [Ys(loc)-100 Ys(loc)+100];
    f1x.YLim = [Xs(loc)-100 Xs(loc)+100];
    set(gcf1, 'Position', [0 800 600 600]);
    % Show same location on mid slice
    f2x.XLim = [Ys(loc)-100 Ys(loc)+100];
    f2x.YLim = [Xs(loc)-100 Xs(loc)+100];
    set(gcf2, 'Position', [800 800 600 600]);
    
    % Draw this nuc on mid slice
    set(0,'CurrentFigure',gcf2);
     Xs = nan;
     Ys = nan;
     for j=1:25
        title({'pick points!!  ', num2str(j)});
        [Xt , Yt] = ginput(1);
        if isempty(Xt)
            break;
        else
            Ys(j) = round(Xt);
            Xs(j) = round(Yt);
            plot(Ys(j), Xs(j),'o');
        end
    end
    Xs = Xs';
    Ys = Ys';
    
    midNucs{i} = [Xs, Ys];

    % write num of nuc in next slice:
    nextNucs(i) = input('nuc num?');
end

save(['matching' num2str(annotSlice)], 'nextNucs');
save(['slice' num2str(annotSlice+5) '_nucs1'], 'midNucs');


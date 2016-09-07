% Main propagate reads the labeled annotated data from one slice,
% First converts it into all nuclei's contours (or snakes)
% Then propagates each snake to the next slices (the n-1 slice and the n+1 slice)
% and runs DP on the surroundings

clc;    % Clear the command window.
clear; % Clear workspace variables.
close all;  % Close all figures

% probabilityImgName = 'probability_maps.tif';
% splitClassesInProbabilityMaps(probabilityImgName);

% Load label (/mask/annotation) image
labelImg = imread('z=1214_label.tif');
sizeImg = size(labelImg);
% Define how dense is the contour sampling:
contourSample = 4;
% lambda defines the normalization in DP minimazation between curvature and and pixal intensity
lambda = 0.9;
slice = 1;

%% Turn label into nuclei contours
nucsContours = label2contour(labelImg, contourSample);


%% For each slice find nucs:

nucsN = size(nucsContours, 1);
% sliceNucs = cell(nucsN ,1);
% 
% save('slice9_Nucs', 'sliceNucs');

load('slice9_Nucs');

screenSize = get(0,'Screensize');

for sliceNum=8:-1:4
    
    % Choose image file
    img = imread('probability_class1.tif', sliceNum);
    % Choose image file
    orgImg = imread('1205_to_1222.tif', sliceNum);
    
    %% Image Filters
    % figure; imshow(img);
    img = imgaussfilt(img, 1);
    thrSobel = 70;
    % Sobel Operator on img with threshold
    [magnitudeImg, directionImg] = imgradient(img, 'sobel');
    magnitudeImg(magnitudeImg < thrSobel) = thrSobel;
    img = -mat2gray(magnitudeImg);
%     figure; imshow(img,[]);
    
    %% run viterbi 
    % to find best match to the nucs in this slice
    % starting location in previous slice
    for i = 10:nucsN
        
        if (~isempty(nucsContours{i}))
            [Xs, Ys, XsInner, XsOuter, YsInner, YsOuter] = addNormThenViterbi(orgImg, img, nucsContours{i}, lambda);
            [Xs, Ys] = resampleSnake4propagate(img, Xs, Ys, contourSample);
        end
        sliceNucs{i} = [Xs, Ys];
        % For testing - Where are the possible bounderies of the nucs
        InnerOuter{i} = [XsInner, XsOuter, YsInner, YsOuter];
    end
   
    %% Show results For correction:

    figure(sliceNum);imshow(orgImg); hold on;
    
    for i=10:nucsN

        temp=sliceNucs{i};
        if ~isempty(temp)
            Xs = temp(:,1);
            Ys = temp(:,2);
            
            loc = round(size(Xs,1)/4);
            xlim([Ys(loc)-150 Ys(loc)+150]);
            ylim([Xs(loc)-150 Xs(loc)+150]);
            set(gcf, 'Position', screenSize); % Maximize figure
            h = plot(Ys, Xs, '--*');
            
            % Manually fix results
            title('Needs fixy?');
            
            [userClickY, userClickX] = ginput(1);
            
            if (abs(userClickX-Xs(1)) < 50)
                % Fix the results of this nucleus
                
                set(h,'Visible','off')
                
                for j=1:25
                    title({'Then fixy!!  ', num2str(j)});
                    [Xt , Yt] = ginput(1);
                    if isempty(Xt)
                        break;
                    else
                        Xs(j) = Xt;
                        Ys(j) = Yt;
                        plot(Xs(j), Ys(j),'o');
                    end
                end
                
                sliceNucs{i} = [Xs, Ys];
            end
        end
        
%          close all;
    end
    
    %% Show results:
%   figure(sliceNum), imshow(orgImg); hold on; impixelinfo;
    for i=10:nucsN
        temp=sliceNucs{i};
        if ~isempty(temp)
            Xs = temp(:,1);
            Ys = temp(:,2);
            h = plot(Ys, Xs, '--*');
        end
    end
    
    nucsContours = sliceNucs;
    save(['slice' num2str(sliceNum) '_Nucs'], 'sliceNucs');
    pause();
    
end


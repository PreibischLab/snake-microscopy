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

%% Turn label into nuclei contours
% nucsContours = label2contour(labelImg, contourSample);
% save('slice9_Nucs', 'nucsContours');

%% For each slice find nucs:
load('slice8_Nucs');
nucsContours = sliceNucs;
nucsN = size(nucsContours, 1);

sliceNucs = cell(nucsN ,6);


for sliceNum=7:-1:1
    
    % Choose image file
     img = imread('probability_class1.tif', sliceNum);
    % Choose image file
    orgImg = imread('1205_to_1222.tif', sliceNum);
    
    %% Image Filters
%     % figure; imshow(img);
%     img1 = imgaussfilt(orgImg, 1);
%     % Sobel Operator on img with threshold
%     [magnitudeImg1, directionImg] = imgradient(img1, 'sobel');
%     magnitudeImg1(magnitudeImg1<80) = 80;
%     magnitudeImg1(magnitudeImg1>350) = 80;
%     img1 = -mat2gray(magnitudeImg1);
%     
%     orgImg2 = orgImg;
%     orgImg2(orgImg2>100) = 0;
%     img2 = imgaussfilt(orgImg2, 1.2);
%     % Sobel Operator on img with threshold
%     [magnitudeImg2, directionImg] = imgradient(img2, 'sobel');
%     magnitudeImg2(magnitudeImg2<40) = 40;
%     magnitudeImg2(magnitudeImg2>350) = 40;
%     img2 = -mat2gray(magnitudeImg2);

    img = imgaussfilt(img, 1);
    % Sobel Operator on img with threshold
    [magnitudeImg, directionImg] = imgradient(img, 'sobel');
    magnitudeImg(magnitudeImg<50) = 50;
    img = -mat2gray(magnitudeImg);
%     figure; imshow(img,[]);


    %% run viterbi 
    % to find best match to the nucs in this slice
    % starting location in previous slice
    
    for i = 1:nucsN    
        % lambda defines the normalization in DP minimazation between curvature and and pixal intensity
        lambda = 0.95;
        
         %% We repeat it 6 times to get optoins for lambda and pick the best segmentation
        for j=1:6
%             img=img1;
        
%             if (~isempty(nucsContours{i,1}))    
            if (size(nucsContours{i,1},1) > 4)
                [Xs, Ys, XsInner, XsOuter, YsInner, YsOuter] = addNormThenViterbi(orgImg, img, nucsContours{i,1}, lambda);
                [Xs, Ys] = resampleSnake4propagate(img, Xs, Ys, contourSample);
            end
            sliceNucs{i, j} = [Xs, Ys];
    %         % For testing - Where are the possible bounderies of the nucs
             InnerOuter{i} = [XsInner, XsOuter, YsInner, YsOuter];
            
             lambda = lambda - 0.05;
        end    
    end

    %% Show results For correction:

    figure(sliceNum);imshow(orgImg); hold on;
    
    for i=1:nucsN
        
        temp=sliceNucs{i, 1};
        
%         if ~isempty(temp)
        if (size(temp,1) >4)
        
            for j=1:6
            
                temp=sliceNucs{i, j};

                Xs = temp(:,1);
                Ys = temp(:,2);
                
                loc = round(size(Xs,1)/4);
                xlim([Ys(loc)-100 Ys(loc)+100]);
                ylim([Xs(loc)-100 Xs(loc)+100]);
                set(gcf, 'Position', [0 0 1400 1400]); % Maximize figure
                
                % Show each result (different lambda) in a different color
                switch j
                    case 1
                        lColor = ':r';
                    case 2
                        lColor = ':y';
                    case 3
                        lColor = ':g';
                    case 4
                        lColor = ':b';
                    case 5
                        lColor = ':m';
                    case 6
                        lColor = ':c';
                        
                end
                
                h = plot(Ys, Xs, lColor, 'LineWidth', 1.5);

            end
            
            % Manually decide which one worked
            title('Which one? 1 2 3 4 5 6 7 for r y g b m c or non?');
            
            whichContour = input('which one?');
            
            if (whichContour==7)
            % Fix the results of this nucleus

                set(h,'Visible','off')

                for j=1:25
                    title({'Then fixy!!  ', num2str(j)});
                    [Xt , Yt] = ginput(1);
                    if isempty(Xt)
                        break;
                    else
                        Ys(j) = Xt;
                        Xs(j) = Yt;
                        plot(Ys(j), Xs(j),'o');
                    end
                end
                sliceNucs{i, 1} = [];
                sliceNucs{i, 1} = [Xs, Ys];
            else
                sliceNucs{i, 1} = sliceNucs{i, whichContour};
            end
            
        end
        
%%      THIS IS THE OLD CODE - Full correction by hand.
%         temp=sliceNucs{i, j};
%         if ~isempty(temp)
%             Xs = temp(:,1);
%             Ys = temp(:,2);
% 
%             loc = round(size(Xs,1)/4);
%             xlim([Ys(loc)-150 Ys(loc)+150]);
%             ylim([Xs(loc)-150 Xs(loc)+150]);
%             set(gcf, 'Position', screenSize); % Maximize figure
%             h = plot(Ys, Xs, ':r');
% 
%             % Manually fix results
%             title('Needs fixy?');
% 
%             [userClickY, userClickX] = ginput(1);
% 
%             if (abs(userClickX-Xs(1)) < 50)
%                 % Fix the results of this nucleus
% 
%                 set(h,'Visible','off')
% 
%                 for j=1:25
%                     title({'Then fixy!!  ', num2str(j)});
%                     [Xt , Yt] = ginput(1);
%                     if isempty(Xt)
%                         break;
%                     else
%                         Xs(j) = Xt;
%                         Ys(j) = Yt;
%                         plot(Xs(j), Ys(j),'o');
%                     end
%                 end
% 
%                 sliceNucs{i} = [Xs, Ys];
%             end
%         end
            
        
        
%       close all;
    end
    
    sliceNucs(:, 2:6) = [];
    
    %% Show results:
%   figure(sliceNum), imshow(orgImg); hold on; impixelinfo;
    for i=10:nucsN
        temp=sliceNucs{i,1};
        if ~isempty(temp)
            Xs = temp(:,1);
            Ys = temp(:,2);
            h = plot(Ys, Xs, '--*');
        end
    end
    
    nucsContours = sliceNucs(:,1);
    save(['slice' num2str(sliceNum) '_Nucs'], 'sliceNucs');
    pause();
    
end


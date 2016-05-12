clc;
clear;
close all;

addpath(genpath('../snake-microscopy'));
templateSize = 64;

img = imread('../DWingPNG/template_affine.png');
figure;imshow(img,[]);
hold on;

%% Landmarks by Yann
if ~(exist('landmarks.mat', 'file') == 2)
    landmarks = [440,503
    983,569
    217,677
    1092,679
    1142,681
    930,686
    1088,713
    938,719
    1105,747
    751,763
    1171,770
    250,775
    769,843
    638,961
    304,557
    332,657
    453,650
    575,648
    697,652
    816,662
    361,876
    491,939
    375,768
    499,760
    626,757
    548,534
    658,559
    768,585
    878,613
    986,645
    579,487
    719,501
    854,529
    761,987
    894,987
    1016,958
    1116,882
    846,744
    877,803
    991,774];
    save('landmarks','landmarks');
else
    load('landmarks.mat');
end

Xs = landmarks(:,1);
Ys = landmarks(:,2);

% Number of landmark points
n = size(Xs, 1);

% Ploting our points:
labels = cellstr(num2str([1:n]'));
scatter(Xs, Ys);
text(Xs, Ys, labels,'FontSize',20);
pause();

%% Make the neighbors matrix - which points adj each point
if ~(exist('adjMat.mat', 'file') == 2)
    adjMat = nan(n,3);
    for i= 1:n
        prompt = ['Who are V' num2str(i) 's neighbors?'];
        adjMat(i,:) = input(prompt);
    end
    save('adjMat','adjMat');
else
    load('adjMat.mat');
end

%% Create the distance matrix between each landmark point to its adjs points
% For 1st derivative
% Plot the connections
children = get(gca, 'children');
delete(children(1:end-1));
adjDist = nan(n,3);
for i= 1:n
    for j= 1:3
        if ~isnan(adjMat(i,j))
            adjX = Xs(adjMat(i,j));
            adjY = Ys(adjMat(i,j));
            adjDist(i,j) = sqrt((Xs(i)-adjX)^2 + (Ys(i)-adjY)^2);
            plot([Xs(i), adjX], [Ys(i), adjY], 'LineWidth', 3);
        end
    end
end
save('adjDist','adjDist');

%% Create the 2nd derivative matrix between each landmark point to its adjs points
% For 1st derivative
% Plot the connections
adjCurve = nan(n,3);
for i= 1:n
    
    adj1X = Xs(adjMat(i,1));
    adj1Y = Ys(adjMat(i,1));
    adj2X = Xs(adjMat(i,2));
    adj2Y = Ys(adjMat(i,2));

    adjCurve(i,1) = (adj1X - 2*Xs(i) + adj2X)^2 + (adj1Y - 2*Ys(i) + adj2Y)^2;
    
    if ~isnan(adjMat(i,3))
        adj3X = Xs(adjMat(i,3));
        adj3Y = Ys(adjMat(i,3));
        adjCurve(i,2) = (adj1X - 2*Xs(i) + adj3X)^2 + (adj1Y - 2*Ys(i) + adj3Y)^2;
        adjCurve(i,3) = (adj2X - 2*Xs(i) + adj3X)^2 + (adj2Y - 2*Ys(i) + adj3Y)^2;
    end
end
save('adjCurve','adjCurve');


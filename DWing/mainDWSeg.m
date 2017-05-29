
%% mainSnake
% Script that sets the snake method used and the parameters
% and start the process by calling the function that will initialize the desired snake 

clc;    % Clear the command window.
clear; % Clear workspace variables.
close all;  % Close all figures

addpath(genpath('../CeleEM'));

%% folder and filename paramters

% Choose image file
folder = '/media/badboy/DATA/taf/test/samples/909_registered/';
folderLandmarks = [folder 'points' filesep];
fileNamePre = 'brightfield_affine';
landmarkPost = '.txt';

templateName = 'template_affine.tif';
landmarkName  = ['template_affine' landmarkPost];
adjacencyName = 'template_affine_points_adj.txt';

snakeGif = 'snakeGif.gif';

%% read landmarks
template  = imread(  [folder          templateName ]);
landmarks = dlmread( [folderLandmarks landmarkName],  ' ');
adjTmp    = dlmread( [folderLandmarks adjacencyName], ' ');


%% VIS **** display template and landmark numberd
clf
figure(1)
imagesc(template)
axis equal tight
colormap('gray')
hold on
scatter(landmarks(:,1),landmarks(:,2), 60, [1, 0, 0],'+')
text(landmarks(:,1)+10,...
     landmarks(:,2)-10,...
     cellstr(num2str((1:size(landmarks,1))')),...
     'Color', [0.6 0 0],... 
     'FontWeight', 'bold')
hold off

%% VIS ****  check graphically the adjacency
for i=1:size(adjTmp,1)
    fprintf('Point %d\n', i);
    clf
    imagesc(template)
    axis equal tight
    colormap('gray')
    hold on
    scatter(landmarks(:,1),landmarks(:,2), 60, [1, 0, 0],'+')
    scatter(landmarks(i,1),landmarks(i,2), 60, [0, 0, 1],'+')
    text(landmarks(:,1)+10,...
         landmarks(:,2)-10,...
         cellstr(num2str((1:size(landmarks,1))')),...
         'Color', [0.6 0 0],... 
         'FontWeight', 'bold')
     hold off
     %%
    for j=2:size(adjTmp,2)
        if adjTmp(i,j) ~= 0
            line(landmarks([adjTmp(i,1) adjTmp(i,j) ],1),...
                 landmarks([adjTmp(i,1) adjTmp(i,j) ],2),...
                 'Color', 'b') 
        end
    end
    pause
end

%% check and reorder ajacency matrix
[m, ix] = ismember(1:size(landmarks,1), adjTmp(:,1));
if sum(~m)>0
    fprintf('landmarks: [ ');
    fprintf('%d ',find(~m));
    fprintf('] not found in adjacency file.\n')
end
adjacency = adjTmp(ix,2:end);
adjacency(adjacency ==0 ) = NaN;

crossPoints = sum(adjacency>0,2) == 3;

%% **** VIS: correct landmark position of the template
clf
imagesc(template)
axis equal tight
colormap('gray')
hPoints = cell(size(landmarks,1),1);
for i=1:size(landmarks,1)
    hPoints{i} = impoint(gca, landmarks(i,:));
end
%% rewrite the landmarks
landmarkNew = zeros(size(landmarks));
for i=1:size(landmarks,1)
    landmarkNew(i,:) = hPoints{i}.getPosition();
end


%% snake parameters
sigma  = 3;   % Gaussian filter sigma (Standard Deviation)

% [alphas betas gammas]
% Alpha - weight of elasticity for list of points (distance)
% Beta  - weight of curvature (abs angular difference)
% Gamma - weight of ext force
semiLdmkParams = [ 0.5   1     0.8 ];
LdmkParams     = [ 0.16   0.3   1   ];

paramPerLdmk = zeros(size(landmarks,1),3);
paramPerLdmk(crossPoints ,:) = repmat(LdmkParams,    sum( crossPoints),1);
paramPerLdmk(~crossPoints,:) = repmat(semiLdmkParams,sum(~crossPoints),1);
paramPerLdmk(14,:) = [0.1, 0.3, 1.2];
paramPerLdmk(1,:) = [0.1, 0.3, 1.2];
templateSize = 64; % template size for template matching

% Choosing hood size to which a point can move to (1->3x3, 2->5x5...)
hoodSize = 50;


%% initialize graphs
figure(1)
% plot template
subplot(1,3,1);
imagesc(template); 
colormap('gray')
set(gca,'xtick',[])
set(gca,'ytick',[])
hold on;
scatter(landmarks(:,1), landmarks(:,2), '+', 'b');
hold off
pos = get(gca, 'Position');
pos(1) = 0;
pos(3) = 0.33;
set(gca, 'Position', pos);
axis equal tight

% plot image with initial landmark position
subplot(1,3,2);
imagesc(img); 
colormap('gray')
set(gca,'xtick',[])
set(gca,'ytick',[])
hold on
scatter(landmarks(:,1), landmarks(:,2), '+', 'r');
hold off

pos = get(gca, 'Position');
pos(1) = 0.33;
pos(3) = 0.33;
set(gca, 'Position', pos)
axis equal tight

% plot image
subplot(1,3,3);
imagesc(img); 
colormap('gray')
set(gca,'xtick',[])
set(gca,'ytick',[])

pos = get(gca, 'Position');
pos(1) = 0.66;
pos(3) = 0.33;
set(gca, 'Position', pos)
axis equal tight

%% extract the smoothed neighborhood of each landmark

templates = createTemplateImages(template, landmarks, 3,32);

%% list images to process
imgList = dir([folder '*.tif']);
imgList = {imgList.name}';
imgList = imgList(~cellfun(@isempty, strfind(imgList,fileNamePre)));

%% Initialize figures (run if needed)
for i=2:5
    figure(i)
    clf
    colormap(map)
end
figure(6)
colormap('gray')
figure(7)
colormap('gray')

%% run the snake %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('+ Processing image             ')
for i=1:numel(imgList)
    str=[num2str(i) '/' num2str(numel(imgList)) ];
    fprintf([repmat('\b',[1 , length(str)+1 ]) str '\n'])
    pause(0.001)
    
    img       = imread(  [folder imgList{i}]);
    ldmkMoving = startIterations(img, landmarks, adjacency, paramPerLdmk,...
                                hoodSize, templates);
                            
    % write the landmark file
    [~,n] = fileparts(imgList{i});
    dlmwrite([folderLandmarks n landmarkPost] ,ldmkMoving,' ');
end 
fprintf('\b -> done\n');
    
%% write a representation of images+ landmark marked in red
M = ...
[ 0 0 1 1 1 1 1 0 0 ;....
  0 0 0 1 1 1 0 0 0 ;....
  0 0 0 1 1 1 0 0 0 ;....
  1 1 1 0 0 0 1 1 1 ;....
  1 1 1 0 0 0 1 1 1 ;....
  1 1 1 0 0 0 1 1 1 ;....
  0 0 0 1 1 1 0 0 0 ;....
  0 0 0 1 1 1 0 0 0 ;....
  0 0 0 1 1 1 0 0 0 ]>0;

for i=1:numel(imgList)
    str=[num2str(i) '/' num2str(numel(imgList)) ];
    fprintf([repmat('\b',[1 , length(str)+1 ]) str '\n'])
    pause(0.001)
    
    [~,n] = fileparts(imgList{i});
    ldmkMoving = dlmread( [folderLandmarks n landmarkPost],  ',');
    img       = imread(  [folder imgList{i}]);
    img = uint8(img);
    I = repmat(img, 1, 1, 3);
    for p=1:size(ldmkMoving,1)
        r = ldmkMoving(p,2)-4:ldmkMoving(p,2)+4;
        c = ldmkMoving(p,1)-4:ldmkMoving(p,1)+4;
        Itmp = I(r,c,1);
        Itmp(M) = 255;
        I(r,c,1) = Itmp;
        Itmp(M) = 0;
        I(r,c,2) = Itmp;
        I(r,c,3) = Itmp;
    end
    imwrite(I, [folder 'test/' imgList{i}])
end

%% rerun the snake on a single image
for i=2:5
    figure(i)
    clf
    colormap(map)
end
figure(6), clf
colormap('gray')
figure(7), clf
colormap('gray')

imgName = 'brightfield_affine00009.tif';
img       = imread(  [folder imgName]);
ldmkMoving = startIterations(img, landmarks, adjacency, paramPerLdmk,...
                            hoodSize, templates);
    
%%
clf
imagesc(img)
axis equal tight
colormap('gray')
hPoints = cell(size(ldmkMoving,1),1);
for i=1:size(ldmkMoving,1)
    hPoints{i} = impoint(gca, ldmkMoving(i,:));
end
%%
figure(1)
clf
imagesc(img)
axis equal tight
colormap('gray')
hold on
scatter(ldmkMoving(:,1),ldmkMoving(:,2), 60, [1, 0, 0],'+')
text(ldmkMoving(:,1)+10,...
     ldmkMoving(:,2)-10,...
     cellstr(num2str((1:size(ldmkMoving,1))')),...
     'Color', [0.6 0 0],... 
     'FontWeight', 'bold')
hold off
    


%% **** VIS: correct landmark position
clf
imgName = 'brightfield_affine00001.tif';
[~,n] = fileparts(imgName);
I = imread([folder imgName]);
ldmkMoving = dlmread( [folderLandmarks n landmarkPost],  ',');
%%
imagesc(img)
axis equal tight
colormap('gray')
hPoints = cell(size(ldmkMoving,1),1);
for i=1:size(landmarks,1)
    hPoints{i} = impoint(gca, ldmkMoving(i,:));
end
%% retrive the landmarks
ldmkMovingNew = zeros(size(ldmkMoving));
for i=1:size(landmarks,1)
    ldmkMovingNew(i,:) = round(hPoints{i}.getPosition());
end
idxNotChanged = find(all(ldmkMovingNew == ldmkMoving,2));

%% actualize the position keeping selected still
ldmkMoving = startIterations(img, landmarks, adjacency, paramPerLdmk,...
                             hoodSize, templates,...
                             idxNotChanged, ldmkMovingNew);

%% rewrite the landmarks

dlmwrite([folderLandmarks n landmarkPost] ,ldmkMoving,' ');

%% gui to modifiy all landmark position

Istack  = cell(numel(imgList),1);
LdmkStack = cell(numel(imgList),2);
for i=1:numel(Istack)
    [~,n] = fileparts(imgList{i});
    LdmkStack{i,1} = dlmread( [folderLandmarks n landmarkPost],  ',');
    LdmkStack{i,2} = LdmkStack{i,1};
    Istack{i} = imread(  [folder imgList{i}]);
end
    
% template

iLdmk = 1;
iImg = 1;

f = figure(1);
clf(f)
hSubs(1) = subplot(1,2,1);
axes(hSubs(1));
colormap('gray')
set(gca,'xtick',[])
set(gca,'ytick',[])

pos = get(gca, 'Position');
pos(1) = 0;
pos(3) = 0.5;
set(gca, 'Position', pos);
axis equal tight

hSubs(2) = subplot(1,2,2);
axes(hSubs(2));
colormap('gray')
set(gca,'xtick',[])
set(gca,'ytick',[])

pos = get(gca, 'Position');
pos(1) = 0.5;
pos(3) = 0.5;
set(gca, 'Position', pos);
axis equal tight

set(f, 'KeyPressFcn', ...
    @(fig_obj , eventDat) ...
    fLdmkChangefigureCallback(fig_obj, eventDat,...
                              hSubs, iLdmk, iImg, 100, ...
                              Istack, LdmkStack, ...
                              template, landmarks));
%% clear

%% rerun sanke with new landmarks

checkedPoints = sum(adjacency>0,2) == 3 || sum(adjacency>0,2) == 1;
fprintf('+ Processing image             ')
for i=1:numel(imgList)
    str=[num2str(i) '/' num2str(numel(imgList)) ];
    fprintf([repmat('\b',[1 , length(str)+1 ]) str '\n'])
    pause(0.001)
%     idxNotChanged = find(all(LdmkStack{i,1} == LdmkStack{i,2},2));
    idxNotChanged = find(checkedPoints);
    
    img       = imread(  [folder imgList{i}]);
    
    ldmkMoving = startIterations(img, landmarks, adjacency, paramPerLdmk,...
                                 hoodSize, templates,...
                                 idxNotChanged, LdmkStack{i,2});
                             
    [~,n] = fileparts(imgList{i});
    dlmwrite([folderLandmarks n landmarkPost] ,ldmkMoving,' ');
end
fprintf('\b-> done\n');
    
    
    
    
    
    
    
    

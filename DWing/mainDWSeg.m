
%% mainSnake
% Script that sets the snake method used and the parameters
% and start the process by calling the function that will initialize the desired snake 

clear; % Clear workspace variables.
close all;  % Close all figures

% addpath(genpath('../CeleEM'));

%% folder and filename paramters

% Choose image file
folder = '/media/badboy/DATA/taf/test/samples/909_registered/';
folderLandmarks = [folder 'points' filesep];
fileNamePre = 'brightfield_affine';
landmarkPost = '.txt';

templateName = 'template_affine.tif';
landmarkName  = ['template_affine' landmarkPost];
adjacencyName = 'template_affine_points_adj.txt';
segmentsName = 'template_affine_points_segments.txt';

snakeGif = 'snakeGif.gif';

%% read landmarks
template  = imread(  [folder          templateName ]);
landmarks = dlmread( [folderLandmarks landmarkName],  ' ');
adjTmp    = dlmread( [folderLandmarks adjacencyName], ' ');
sgmTmp    = dlmread( [folderLandmarks segmentsName], ' ');
nLdmk = size(landmarks,1);
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

%% check and reorder adjacency matrix
[m, ix] = ismember(1:size(landmarks,1), adjTmp(:,1));
if sum(~m)>0
    fprintf('landmarks: [ ');
    fprintf('%d ',find(~m));
    fprintf('] not found in adjacency file.\n')
end

adjacency = adjTmp(ix,2:end);
adjacency = arrayfun(@(x) adjacency(x,adjacency(x,:)~=0), (1:size(adjacency,1))','Unif',false); 

ldmkClass = zeros(size(adjacency,1),1);
tmp = cellfun(@numel, adjacency);
crossPoints = tmp == 3;
ldmkClass(crossPoints) = 1;
ldmkClass(tmp == 1) = 2; %endpoints
ldmkClass(tmp == 2) = 3;% semi

%% check and reorder segments matrix

segments = arrayfun(@(x) sgmTmp(x,sgmTmp(x,:)~=0), (1:size(sgmTmp,1))','Unif',false); 
% set the points equidistant
segmentsPos = cell(size(segments));
for i=1:numel(segments)
    segmentsPos{i} = linspace(0,1,numel(segments{i}));
end

%% get neighborhood for landmarks from semgents
bigAdjTmp = cell2mat(cellfun(@(x) [x(1) x(end)],segments,'Unif',false));
LbigDist = cell(size(landmarks,1),1);
LbigAngle = cell(size(landmarks,1),1);
for i=1:size(landmarks,1)
    LbigDist{i}  = setdiff(unique(bigAdjTmp(sum(bigAdjTmp==i,2)>0,:)),i);
    if numel(LbigDist{i})>1
        LbigAngle{i} = [LbigDist{i}(1:end-1) LbigDist{i}(2:end)];
    end
end

%% **** GUI: correct landmark position of the template
clf
imagesc(template)
axis equal tight
colormap('gray')
hPoints = cell(size(landmarks,1),1);
for i=1:size(landmarks,1)
    hPoints{i} = impoint(gca, landmarks(i,:));
end
%% **** GUI: rewrite the landmarks
landmarkNew = zeros(size(landmarks));
for i=1:size(landmarks,1)
    landmarkNew(i,:) = hPoints{i}.getPosition();
end

%% base paramters

templateSize = 60; % template size for template matching
sigma  = 3;   % Gaussian filter sigma (Standard Deviation)

% Choosing hood size to which a point can move per iteration
% (1->3x3, 2->5x5... 30-> 61x61)
hoodSize = 50;
bLearnFromStack = true(1);

%% extract the smoothed neighborhood of each landmark
if bLearnFromStack
    fol = [folder 'wrapped/'];
    listWrapped = dir([fol '*.tif']);
    listWrapped = {listWrapped.name}';
    Istack = cell(1,1,numel(listWrapped));
    G = fspecial('gaussian',[5 5],sigma);
    for i=1:numel(listWrapped)
        Istack{i} = imread([fol listWrapped{i}]);
        Istack{i} = imfilter(Istack{i},G,'same');
    end
    Istack = cell2mat(Istack);
    templates = createTemplateImages(mean(Istack,3), landmarks, 0.01, templateSize);
%     figure
%     imagesc(mean(Istack,3))
else
    templates = createTemplateImages(template, landmarks, sigma, templateSize);
end

%%
fLearnTemplate(Istack, landmarks, ldmkClass, sigma);

%% read landmarks stack
if bLearnFromStack
    fol = [folderLandmarks 'corrected/'];
    txtList = dir([fol 'brightfield*.txt']);
    txtList = {txtList.name}';
    landmarksStack = cell(numel(txtList),1);
    for i=1:numel(txtList)
        landmarksStack{i} = dlmread( [fol txtList{i}],  ',');
    end
    landmarksStack = reshape(landmarksStack,1,1,[]);
    landmarksStack = cat(3,landmarksStack{:});
else
    % take the template
    landmarksStack = landmarks;
end
%% compute the shape model based on the template landmarks
% [alphas betas gammas]
% Alpha - weight of elasticity for list of points (distance)
% Beta  - weight of curvature (abs angular difference)
% Gamma - weight of ext force
weights = [ 0.2   0.2   1   ;...% LdmkParams
            0.5   1     0.8 ;...% end points
            1     0.5   1.5  ]; % semiLdmkParams
if bLearnFromStack
    shapeModel = fEstimateShapeModel(landmarksStack, ldmkClass, weights,...
                                     segments, segmentsPos, adjacency);
    
else
    shapeModel = fEstimateShapeModel(landmarksStack, ldmkClass, weights,...
                                     segments, segmentsPos,...
                                     adjacency, LbigDist, LbigAngle);
end
%% VIS***:show the neighbor selection
for i=1:size(landmarks,1)
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
     adjC = shapeModel.curvIdx{i};
     adjD = shapeModel.distIdx{i};
     col = {'b','r','g','m'};
%     for j=1:size(adjC,1)
%         line(landmarks([i adjC(j,1) ],1)+4*(j-1),...
%              landmarks([i adjC(j,1) ],2)+4*(j-1),...
%              'Color', col{j}) 
%         line(landmarks([i adjC(j,2) ],1)+4*(j-1),...
%              landmarks([i adjC(j,2) ],2)+4*(j-1),...
%              'Color', col{j}) 
%     end
    for j=1:numel(adjD)
        line(landmarks([i adjD(j) ],1),...
             landmarks([i adjD(j) ],2),...
             'Color', col{j})
    end
    pause
end


%% list images to process
imgList = dir([folder '*.tif']);
imgList = {imgList.name}';
imgList = imgList(~cellfun(@isempty, strfind(imgList,fileNamePre)));

%% Initialize figures (run if needed)
figure(2)
clf
for i=1:4
    h = subplot(2,3,i);
    colormap(h,map)
end
h = subplot(2,3,5);
colormap(h, 'gray')
h = subplot(2,3,6);
colormap(h, 'gray')

%% run the snake %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folTmp = '/media/badboy/DATA/taf/test/samples/909_registered/test_vis/';
fprintf('+ Processing image             ')
for i=1:numel(imgList)
    str=[num2str(i) '/' num2str(numel(imgList)) ];
    fprintf([repmat('\b',[1 , length(str)+1 ]) str '\n'])
    pause(0.001)
    
    img       = imread(  [folder imgList{i}]);
    [ldmkMoving, ~] = startIterations(img, shapeModel,...
                                        templates , hoodSize  );
                            
    % write the landmark file
    [~,n] = fileparts(imgList{i});
    dlmwrite([folderLandmarks n landmarkPost] ,ldmkMoving,' ');
end 
fprintf('\b -> done\n');
%% run the for the output  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('+ Processing image             ')
for i=1:numel(imgList)
    str=[num2str(i) '/' num2str(numel(imgList)) ];
    fprintf([repmat('\b',[1 , length(str)+1 ]) str '\n'])
    pause(0.001)
    
    img       = imread(  [folder imgList{i}]);
    [ldmkMoving, out] = startIterations(img, shapeModel, landmarks,...
                                        templates , hoodSize    );
                            
    % write the landmark file
    [~,n] = fileparts(imgList{i});
    f = fieldnames(out);
    for j=1:numel(f)
        Iaff = out.(f{j});
        
        imwrite(Iaff,mapTmp,[folderLandmarks n landmarkPost])
        
    end
    dlmwrite([folderLandmarks n landmarkPost] ,ldmkMoving,' ');
end 
fprintf('\b -> done\n');
%% write a representation of images + landmark marked in red
M = ...
[ 0 0 1 1 1 1 1 0 0 ;....
  0 0 0 1 1 1 0 0 0 ;....
  1 0 0 1 1 1 0 0 1 ;....
  1 1 1 0 0 0 1 1 1 ;....
  1 1 1 0 0 0 1 1 1 ;....
  1 1 1 0 0 0 1 1 1 ;....
  1 0 0 1 1 1 0 0 1 ;....
  0 0 0 1 1 1 0 0 0 ;....
  0 0 1 1 1 1 1 0 0 ]>0;

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

%% Procruste superimposition
% reads the landmarks from file
LdmkStack = cell(numel(imgList),1);
for i=1:numel(LdmkStack)
    [~,n] = fileparts(imgList{i});
    LdmkStack{i} = dlmread( [folderLandmarks n landmarkPost],  ',');
end

% run the procruste superimposition
[SupDissim,SupData, SupTform] =  fProcrustesSupp(LdmkStack, false);

% average shape


%% point pairs
PPairs = [reshape(repmat((1:size(adjacency))',1,size(adjacency,2)),[],1) reshape(adjacency,[],1)];
PPairs = PPairs(~any(isnan(PPairs),2),:);
PPairs = sort(PPairs,2);
PPairs = unique(PPairs, 'rows');

%% plot shape for all individuals (procruste)
clf
mapD = fDistingColors(40,[1 1 1]);
clf
for i=1:numel(SupData)
    for j=1:size(PPairs,1)
        line(SupData{i}(PPairs(j,:),1),SupData{i}(PPairs(j,:),2),'Color',mapD(i,:))
    end
end
set(gca,'Ydir','reverse')
axis equal tight
hl = axis;
hl([1 3]) = hl([1 3]) - d;
hl([2 4]) = hl([2 4]) + d;
axis(hl)
%% plot points + model  (procruste)
Mldmk = cell2mat(reshape(SupData,1,1,[]));
avLdmks = mean(Mldmk,3);
clf
for i=1:numel(SupData)
    scatter(SupData{i}(:,1),SupData{i}(:,2),50, mapD(i,:), '+')
    hold on
end
d = diff(reshape(hl,2,2),1, 1)*0.1;
set(gca,'Ydir','reverse')
axis equal tight
hl = axis;
hl([1 3]) = hl([1 3]) -d;
hl([2 4]) = hl([2 4]) +d;
axis(hl)

% draw adjacency lines
for i=1:size(PPairs,1)
    line(avLdmks(PPairs(i,:),1),avLdmks(PPairs(i,:),2), 'Color', [0 0 0])
end
%% plot points + model for initial registration (initReg)
clf
Mldmk = cell2mat(reshape(LdmkStack,1,1,[]));
avLdmks = mean(Mldmk,3);

for i=1:numel(LdmkStack)
    scatter(LdmkStack{i}(:,1),LdmkStack{i}(:,2),50, mapD(i,:), '+')
    hold on
end
d = diff(reshape(hl,2,2),1, 1)*0.1;
set(gca,'Ydir','reverse')
axis equal tight
hl = axis;
hl([1 3]) = hl([1 3]) -d;
hl([2 4]) = hl([2 4]) +d;
axis(hl)

% draw adjacency lines
for i=1:size(PPairs,1)
    line(avLdmks(PPairs(i,:),1),avLdmks(PPairs(i,:),2), 'Color', [0 0 0])
end
%% plot shape for all individuals (initReg)
clf
mapD = fDistingColors(40,[1 1 1]);
clf
for i=1:numel(LdmkStack)
    for j=1:size(PPairs,1)
        line(LdmkStack{i}(PPairs(j,:),1),LdmkStack{i}(PPairs(j,:),2),'Color',mapD(i,:))
    end
end
set(gca,'Ydir','reverse')
axis equal tight
hl = axis;
hl([1 3]) = hl([1 3]) -d;
hl([2 4]) = hl([2 4]) +d;
axis(hl)

%% rerun the snake on a single image
figure(2)
clf
for i=1:4
    h = subplot(2,3,i,'replace');
    colormap(h,map)
end
h = subplot(2,3,5);
colormap(h, 'gray')
h = subplot(2,3,6);
colormap(h, 'gray')
%%
imgName = 'brightfield_affine00004.tif';
img       = imread(  [folder imgName]);
[ldmkMoving, out, ldmkHist] = startIterations( img, shapeModel,...
                                     templates, hoodSize   );
%%
figure(1)
colormap(map)
imagesc(out.total), axis equal tight
hold on
scatter(ldmkMoving(1:end-1,1),ldmkMoving(1:end-1,2),70, 'r','+')
scatter(ldmkMovingOld(1:end-1,1),ldmkMovingOld(1:end-1,2),70, 'w','+')
hold off
%%

figure(1)
colormap('gray')
imagesc(img), axis equal tight
hold on
for i=1:size(ldmkHist,1)
    plot(reshape(ldmkHist(i,1,:),[],1),...
         reshape(ldmkHist(i,2,:),[],1),...
         'r','LineWidth',2)
end
hold off

%% draggalbe points on the last created landmarks
clf
imagesc(img)
axis equal tight
colormap('gray')
hPoints = cell(size(ldmkMoving,1),1);
for i=1:size(ldmkMoving,1)
    hPoints{i} = impoint(gca, ldmkMoving(i,:));
end

%% ***** VIS: plot the numbered landmarks
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
imgName = 'brightfield_affine00004.tif';
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
ldmkMoving = startIterations(img, shapeModel, ldmkMovingNew,...
                                   templates, hoodSize,... 
                                  idxNotChanged);

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
    
    ldmkMoving = startIterations(img, shapeModel, LdmkStack{i,2},...
                                      templates , hoodSize,... 
                                      idxNotChanged);
                             
    [~,n] = fileparts(imgList{i});
    dlmwrite([folderLandmarks n landmarkPost] ,ldmkMoving,' ');
end
fprintf('\b-> done\n');
    
    
    
    
    
    
    
    

%% Snake Implementations
%  KASS and GREEDY snakes

clc;    % Clear the command window.
clear; % Clear workspace variables.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.

%% Parameters

% Type of snake implementation - either 'KASS' or 'GREEDY'
method = 'KASS';
% Number of snake points
nPoints = 200;
% Choose num of iterations
iterations = 40000;
% Choosing if starting snake will be an automatic cyrcle or a manual input
% from user - either 'CIRCLE' or 'MANUAL'
StartingSnakeMethod = 'MANUAL';
% Gaussian filter sigma (Standard Deviation)
sigma = 1;
% In Greedy only - Choosing neighborhood size - in which any point in the snake can move to.
% for 1 the neighborhood will be a 3X3, for 2 5X5...
hoodSize = 10;
% Choosing alpha - weight of elasticity 
alpha = 0.01;
% Choosing beta - weight of curvature 
beta = 0.05;
% Only in greedy - Choosing gamma - weight of gradient (sobol)
gamma = 5;
% Choose delta - step size - Only in Kass.
delta = 1;
% Choose threshold on sobel before normalization:
threshSobel = 50;
% Choose image file
img = imread('n2copy.jpg');

% Getting img dimentions
[rows, columns] = size(img);

%% Create the initial snake

% If user chose to manually select the initial snake
if strcmp(StartingSnakeMethod, 'MANUAL')
    % Number of points user has to insert - those points will be interpulated 
    n=10;
    snake = setManualSnake(img, n, nPoints);
    % Since after the interpolation snake size might not be exactly nPoints (but very close):
    nPoints = size(snake,2);
% Else - if initial snake is set by a circle     
else
    % For initial snake contour - making a circle
    radius = floor(rows/2) - 10;
    center = [radius+10, radius+10];
    
    snake = zeros(2,nPoints);

    for i= 1:nPoints
        snake(1,i) = center(1) + floor(radius*cos(i*2*pi/nPoints) + 0.5);
        snake(2,i) = center(2) + floor(radius*sin(i*2*pi/nPoints) + 0.5);
    end
    
    figure;
    imshow(img);
    hold on;
end

%% Pre-processing 

% Show initial snake
plot(snake(1,:),snake(2,:),'*');
impixelinfo;

orgImg = img;
% Gausian Filter on img
img = imgaussfilt(img, sigma);

% Sobel Operator on img with threshold
[magnitudeImg, directionImg] = imgradient(img, 'Sobel');
magnitudeImg(magnitudeImg<threshSobel) = threshSobel;
img = -mat2gray(magnitudeImg);

% % SOBEL FOR GREEDY - comment last line, uncomment this line
% img = 254*mat2gray(magnitudeImg);

% Initializing snake parameter
snake(3,:) = alpha;
snake(4,:) = beta;

%% Call snake fuction - Kass or Greedy 
if strcmp(method, 'GREEDY')
    % Making image frame - against out of bounderies
    newImg = zeros(rows+100, columns+100);
    newImg(51:end-50,51:end-50) = img;
    snake(1,:) = snake(1,:)+50;
    snake(2,:) = snake(2,:)+50;
    
    % Define weight for external forces 
    snake(5,:) = gamma;
   
    snake = greedy(newImg, snake, iterations, hoodSize, orgImg);
else
    snake = kass(img, snake, iterations, delta, orgImg, sigma, StartingSnakeMethod);
end
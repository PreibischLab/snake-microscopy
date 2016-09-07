%% mainSnake
% Script that sets the snake method used and the parameters
% and start the process by calling the function that will initialize the desired snake 

clc;    % Clear the command window.
clear; % Clear workspace variables.
close all;  % Close all figures

%% Setting the snakes parameters
%slice = 254;
% Choose image file
%img = imread(['/Users/ebahry/Desktop/images_em_head/z' num2str(slice) '.tif'],'PixelRegion',{[1250 5000],[3500 8800]});
% if (ndims(img)>2)
%     img = rgb2gray(img);
% end
% img = imresize(img,0.4);
% [img, rect] = imcrop(img);
% imwrite(img, 'cell.tiff');
%img=img(3200:8800,1000:5000);
% imshow(img);impixelinfo;

% Choose image file
img = imread(['../classification_result.tif'],1);

% Type of snake implementation - KASS, BALLOON, DUAL
method = 'BALLOON';
% Initial snake shape - either 'MANUAL' or radius of circle (e.g. 50) 
% If dual - array with 2 values - first kass and then for balloon snake
initShape = {'50'};
% Gaussian filter sigma (Standard Deviation)
sigma = 3.5;
% Threshold on sobel (before normalization):
thrSobel = 30;
% Alpha - weight of elasticity:
alpha = [0.015];
% Beta - weight of curvature 
beta = [0.01];

% Rho - weight of expention (for Balloon & dual)
rho = 0.2;

% Delta - step size (not used in greedy).
delta = 1;

% In Greedy only:
% Choosing hood size to which a point can move to (1->3x3, 2->5x5...)
hoodSize = 10;
% Gamma - weight of gradient (sobol)
gamma = 5;

% NEED TO ADD
% multiSnakeMode = 'NO'; % some of the function already exists
%(setMultiInitSnakes(img, nPoints));

%% Calling the initialization of the snake
initSnake(img, method, initShape, sigma, thrSobel, alpha, beta, rho, delta, hoodSize, gamma, slice);
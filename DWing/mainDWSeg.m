%% mainSnake
% Script that sets the snake method used and the parameters
% and start the process by calling the function that will initialize the desired snake 

clc;    % Clear the command window.
clear; % Clear workspace variables.
close all;  % Close all figures

addpath(genpath('../CeleEM'));

%% Setting the snakes parameters

% Choose image file
folder = '../../DWingPNG/';
fileNamePre = 'brightfield_affine0000';
% 
% for i=1:23
%     if i==10
%         fileNamePre = fileNamePre(1:end-1);
%     end
%     img = imread([folder fileNamePre num2str(i) '.tif']);
%     figure; imshow(img,[]);
%     pause();
% end

% img = imread([folder 'template_affine.tif']);
% figure; imshow(img,[]);

img = imread([folder 'brightfield_affine00000.png']);
% NEED TO ADD - IMAGE RESIZE: img = imresize(img, 0.25);
figure; imshow(img,[]); hold on;

% Gaussian filter sigma (Standard Deviation)
sigma = 2.5;
% Alpha - weight of elasticity for list of points:
alphas = [0.001 0.01];
% Beta - weight of curvature 
betas = [0.000055 0.015];
% Gamma - weight of ext force
gammas = 5;

% template size for template matching:
templateSize = 32;
% Choosing hood size to which a point can move to (1->3x3, 2->5x5...)
hoodSize = 10;

%% Calling the initialization of the snake
initSnake(img, sigma, alphas, betas, gammas, templateSize, hoodSize);
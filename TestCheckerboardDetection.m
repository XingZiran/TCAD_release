clear;
clc;
close all;
addpath('./MyDetector/');
addpath('./Pattern/');
%% load pattern
load('PatternInfo.mat');
figure;imshow(Pattern);hold on;
for i = 1 : size(PatternPts,1)
    label = sprintf('%d',i);
    text(PatternPts(i,1), PatternPts(i,2), label,'BackgroundColor', [1 1 1]);
end
title('Pattern');
%% load photos
Image = imread('./data/000973.jpg');
Image = rgb2gray(Image);

% To avoid long processing time. Uncomment these lines when your images are large.
% [r,c,ch] = size(Image);
% wc = 1920;
% if c > wc
%     Image = imresize(Image,wc/c);
% end
% Image = imresize(Image,[1080,1920]);
% Image = imresize(Image, 2);
%% detect image points of pattern : matlab
tic;
Matlab_Pts = detectCheckerboardPoints(Image);
toc;

if ~isempty(Matlab_Pts)
    figure;imshow(Image);hold on;plot(Matlab_Pts(:,1),Matlab_Pts(:,2),'ro');
    for i = 1 : size(Matlab_Pts,1)
        label = sprintf('%d',i);
        text(Matlab_Pts(i,1), Matlab_Pts(i,2), label, 'BackgroundColor', [1 1 1]);
    end
else
    figure;imshow(Image);
    disp('Matlab failed !');
end
title('Matlab result');
%% detect image points of pattern and extract descriptors
tic;
[I_Pts,boardSize] = detectMyPatternPoints(Image,false);
toc;
if ~isempty(I_Pts)
    ID = findPtsID( I_Pts,boardSize,squareSize,Pattern,PatternPts,PatternMatrixSize,Image,false);
    figure;imshow(Image);hold on;
    for i = 1 : size(I_Pts,1)
        label = sprintf('%d',ID(i));
%         label = sprintf('%d',i);
        text(I_Pts(i,1), I_Pts(i,2), label,'BackgroundColor', [1 1 1]);
%         plot(I_Pts(i,1), I_Pts(i,2),'ro');
    end
    disp('Detection complete!');
else
    disp('No checkerboard detected!');
end
title('TCAD result');
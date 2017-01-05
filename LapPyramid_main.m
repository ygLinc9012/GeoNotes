% ================================
% Demo of Image Laplace Pyramid
% This method is based on the following paper:
%
% Burt, Peter, and Edward Adelson. 
% "The Laplacian pyramid as a compact image code." 
% IEEE Transactions on communications 31.4 (1983): 532-540.
% ================================

clc; clear all; close all;
imgOrig = im2double(imread('img.bmp'));
maxLevel = 4;
sigma = 2;

% Reduce
img = imgOrig;
for level = 1:maxLevel
    temp = imgaussfilt(img, sigma);
    idx = [1:2:size(img,1)];
    L{level} = temp(idx,idx);
    img = L{level};
end

% Expand
for level = maxLevel:-1:1
    img = L{level};
    temp = imresize(img,2);
    expand{level} = imgaussfilt(temp, sigma);
end

% Laplace Pyramid
img = imgOrig;
for level = 1:maxLevel-1
    B{level} = img - expand{level};
    img = L{level};
    
    % plot
    figure(2)
    subplot(1,maxLevel,level)
    imshow(B{level}+0.5)
    title(strcat('lap pyr lv', int2str(level)));
end
subplot(1,maxLevel,maxLevel)
imshow(L{maxLevel})
title(strcat('lap pyr lv', int2str(maxLevel)));

% Reconstruction test
figure(3)
subplot(1,2,1)
imshow(imgOrig)
title('original image')
subplot(1,2,2)
imshow(B{1} + expand{1})
title('reconstruct image')


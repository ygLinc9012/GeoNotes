% ================================
% Demo of Texture Optimization using Markov Random Field energy
% This method is based on the section 3 of the following paper:
%
% Kwatra, Vivek, et al. 
% "Texture optimization for example-based synthesis." 
% ACM Transactions on Graphics (ToG) 24.3 (2005): 795-802.
% ================================

% Initialization
clc; clear all; close all;
inputImg = im2double(rgb2gray(imread('texture.png'))); % Input image Z
[inputH, inputW] = size(inputImg);
figure(1)
imshow(inputImg);

% Parameter
maxIter = 9;
outputW = 100; % Output Image width
outputH = 100; % Output Image height
w = 6; % window size
sampleRate = ceil(w/2); % sample rate on output
k = (2*w+1)^2; % number of pixel in a window
NNF = zeros((floor((outputH-2*w-1)/sampleRate)+1) *(floor((outputW-2*w-1)/sampleRate)+1), 1);
outputImg = rand(outputH,outputW); % initialize random noise output image

% Construct input patch data for Nearest Neighbor computation
inputPatch = zeros((inputH-2*w)*(inputW-2*w), k);
for i = 1 : inputH-2*w
    for j = 1 : inputW-2*w
        idx = (i-1)*(inputW-2*w) + j;
        inputPatch(idx, :) = reshape(inputImg(i:i+2*w, j:j+2*w, :), 1, k);
    end
end

% Start Optimization
for iter = 1:maxIter
    
    % Step1: Update NNF
    NNF_old = NNF;
    NNFIdx = 1;
    totalEnergy = 0;
    for i = w+1:sampleRate:outputH-w
        for j = w+1:sampleRate:outputW-w
            outputPatch = reshape(outputImg(i-w:i+w, j-w:j+w, :), 1, k); % patch on X
            dist = sum(abs(inputPatch - repmat(outputPatch,[size(inputPatch,1),1])),2);
            [dist,idx] = min(dist);
            totalEnergy = totalEnergy + dist;
            NNF(NNFIdx) = idx;
            NNFIdx = NNFIdx + 1;
        end
    end
    fprintf('total MRF energy = %f\n', totalEnergy);
    
    % Check terminate condition
    if (NNF_old == NNF)
        break;
    end
    
    % Step2: Update Output Image
    outputNew = zeros(size(outputImg));
    inputNNF = inputPatch(NNF,:);
    XIdx = 1;
    overlapCount = zeros(size(outputImg));
    for i = w+1:sampleRate:outputH-w
        for j = w+1:sampleRate:outputW-w
            overlapCount(i-w:i+w, j-w:j+w, :) = overlapCount(i-w:i+w, j-w:j+w, :) + 1;
            temp = reshape(inputNNF(XIdx,:), [sqrt(k), sqrt(k)]);
            outputNew(i-w:i+w, j-w:j+w, :) = outputNew(i-w:i+w, j-w:j+w, :) + temp;
            XIdx = XIdx + 1;
        end
    end
    
    % Plot output images
    outputNew = outputNew ./ overlapCount;
    figure(2)
    subplot(ceil(sqrt(maxIter)),ceil(sqrt(maxIter)),iter)
    imshow(outputImg)
    title(strcat('iter = ',num2str(iter)));
    
    outputImg = outputNew;
end








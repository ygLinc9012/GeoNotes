% ================================
% Demo of computing mesh saliency using spectral methods
% This method is based on the following paper:
%
% Song, Ran, et al. 
% "Mesh saliency via spectral processing." 
% ACM Transactions on Graphics (TOG) 33.1 (2014): 6.
% ================================

clc; clear all; close all;
numFace = 3000;

% Compute Laplace eVal and eVec
[V,F] = readOBJ('lion.obj') ;
[F,V] = reducepatch(F,V,numFace);
L = CotanLaplace(V, F);

[eVec, eVal] = eig(full(L));
eVal = diag(eVal);
logEVal = log(abs(eVal));

% Compute A
paddEVal = padarray(logEVal,4,'replicate','both');
A = filter(ones(1,9)/9,1,paddEVal);
A(1:4) = [];
A(end-3:end) = [];
R = abs(logEVal - A);
R = diag(exp(R));

% Compute W (euclidean distance mat)
adjMat = triangulation2adjacency(F);
[i,j] = find(adjMat);
dist = sum((V(i,:) - V(j,:)).^2, 2);
W = sparse(i,j,dist);  
W(W>0) = 1./W(W>0);
W = (W+W')/2;

% Compute Saliency
S = eVec * R * eVec' .* full(W);
S = S';
saliency = sum(S,2);

% Plot
trimesh(F,V(:,1),V(:,2),V(:,3),log(abs(saliency)),'FaceColor','interp','edgecolor','interp')
axis off
axis equal
shading interp
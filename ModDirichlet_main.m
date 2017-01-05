% ================================
% Demo of Modified Dirichlet Energy (MDE) Operator
% This demo function computes the MDE operator 
% and visualize its eigenvalue and eigenvectors
% 
% Reference:
% Hildebrandt, Klaus, et al. 
% "Modal shape analysis beyond Laplacian." 
% Computer Aided Geometric Design 29.5 (2012): 204-218.
% ================================

clc; clear all; close all;
[V,F] = readOBJ('fandisk.obj');
MDE = ModDirichlet(V, F);

% Compute Mass Matrix for generalized eigenvalue problem
[~, vertAreaList, ~, ~] = ComAreaAndNormal(F,V);
Mv = sparse(size(vertAreaList,1), size(vertAreaList,1)); % mass matrix
for idx = 1:size(vertAreaList,1)
    %Mv = diag(vertAreaList);
    Mv(idx, idx) = vertAreaList(idx);
end

% Compute eigenvalue and eigenvectors
numEigs = 6;
[eVec, eVal] = eigs(MDE + 1e-6.*speye(size(MDE,1),size(MDE,1)), Mv, numEigs, 'sm');
[~ ,i] = sort(diag(eVal)); % sort
eVal = sum(eVal,2);
eVal = eVal(i) - 1e-6;
eVec = eVec(:,i);

% Visualize
figure(1)
plot(eVal)
title('eigenvalue')
grid on

figure(2)
for i = 1:numEigs
    subplot(2,ceil(numEigs/2),i)
    trimesh(F,V(:,1),V(:,2),V(:,3),eVec(:,i),'FaceColor','interp','edgecolor','interp');
    colormap jet
    axis equal
    axis off
    title(strcat('eigenvector#',int2str(i)))
end


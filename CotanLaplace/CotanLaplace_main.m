% ================================
% Demo of Cotan Laplace Operator
% This demo function compute cotangent weighted laplace operator 
% and visualize its eigenvalue and eigenvectors
%
% Reference:
% http://ddg.cs.columbia.edu/SGP2014/LaplaceBeltrami.pdf
% ================================


clc; clear all; close all;
[V,F] = readOBJ('bunny.obj');
L = CotanLaplace(V, F);
[~, vertAreaList, ~, ~] = ComAreaAndNormal(F,V);
Mv = sparse(size(vertAreaList,1), size(vertAreaList,1)); % mass matrix
for idx = 1:size(vertAreaList,1)
    Mv(idx, idx) = vertAreaList(idx);
end

% Compute eigenvalue and eigenvectors
numEigs = 6;
[eVec, eVal] = eigs(L + 1e-6.*speye(size(L,1),size(L,1)), Mv, numEigs, 'sm');
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
    view(0,90)
    title(strcat('eigenvector#',int2str(i)))
end


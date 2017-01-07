% ================================
% Demo of computing Dirac sphere
% This method is based on the section 5 of the following paper:
%
% Crane, Keenan, Ulrich Pinkall, and Peter Schröder. 
% "Spin transformations of discrete surfaces." 
% ACM Transactions on Graphics (TOG). Vol. 30. No. 4. ACM, 2011.
% ================================

clc; clear all; close all;
[V,F] = readOBJ('Sphere.obj');
D = Dirac(F,V);
nV = size(V,1);
nF = size(F,1);
[faceAreaList, ~, ~, ~] = ComAreaAndNormal(F,V);

%% Compute square of Dirac operator
% face area matrix
colIdx = [1:nF];
colIdx = reshape(colIdx, [], 1);
colIdx = [colIdx*4-3; colIdx*4-2; colIdx*4-1; colIdx*4];
Mf = sparse(colIdx, colIdx, [faceAreaList; faceAreaList; faceAreaList; faceAreaList], 4*nF, 4*nF);
% face-vert matrix
colIdx = [1:nF];
colIdx = [colIdx; colIdx; colIdx];
colIdx = reshape(colIdx, [], 1);
colIdx = [colIdx*4-3; colIdx*4-2; colIdx*4-1; colIdx*4];
rowIdx = reshape(F', [], 1);
rowIdx = [rowIdx*4-3; rowIdx*4-2; rowIdx*4-1; rowIdx*4];
B = sparse(colIdx, rowIdx, 1/3, 4*nF, 4*nV);
% Dirac Square
D2 =  D' * Mf * D;
D2dif = max(max(D2' - D2)); % check hermitian (if this value is small, D2 is hermitian)
fprintf(strcat('Dirac square difference = ', num2str(D2dif)));
D2 = (D2' + D2) ./ 2; % make it hetmitian

%% Build Dirac mass matrix
massMat = B'*Mf*D;
massMat = (massMat' + massMat) ./ 2; % make it real hetmitian

%% Eigs
[eVec, eVal] = eigs(D2+1e-6.*speye(size(D2,1),size(D2,2)), massMat+1e-6.*speye(size(D2,1),size(D2,2)), 100, 'sm');
[~ ,i] = sort(diag(abs(eVal))); % sort
eVal = real(sum(eVal,2));
eVal = eVal(i);
eVec = eVec(:, i);

% Plot eigenvalue
figure(1)
scatter([1:100], eVal) % eigenvalue and eigenvector have duplication because of quaternion space
grid on
xlabel('index')
ylabel('eigen values')

%% Dirac Sphere
i = 1;
color = ones(size(V,1),1);
for idx = [1 5 9 17 21 25 29 35 37]
    subplot(3,3,i)
    Vnew = SpinTrans(F,V,real(eVec(:,idx)));
    trimesh(F,Vnew(:,1),Vnew(:,2),Vnew(:,3), color, 'edgealpha', 0.4,'facealpha',0.4)
    axis equal
    axis off
    shading interp  
    i = i + 1;
end

% ================================
% Demo of Dirac Disk
% This method is based on the section 6.5.2 in the following paper:
%
% Crane, Keenan Michael. 
% Conformal geometry processing. 
% Diss. California Institute of Technology, 2013.
% ================================

clc; clear all; close all;

[V,F] = readOBJ('disk.obj');
nF = size(F, 1); % number of faces
nV = size(V, 1); % number of vertices
BV = compute_boundary(F); % boundary vertex indices
BV = sort(BV');

% Reorder boundary vertices
V_b = V(BV, :);
V(BV,:) = [];
V = [V; V_b];
idxMat = [1:size(V,1)]';
idxMat(BV) = [];
idxMat = [idxMat; BV];
Fnew = F;
for i = 1:size(V,1)
    Fnew(find(F == i)) = find(idxMat == i);
end
F = Fnew;
BV = [size(V,1)-size(BV,1)+1 : size(V,1)];
[faceAreaList, vertAreaList, FN, VN] = ComAreaAndNormal(F,V);

% Build Dirac oprator
D = Dirac(F, V);

% Compute Binormal vectors
% 1. Boundary normals
VN = VN(BV,:); 
% 2. Boundary tangent
adjV = adjacency_list(F,BV);
for idx = 1:size(adjV,1)
    i = BV(idx);
    adjBV = intersect(BV, adjV{idx});
    if size(adjBV,2) ~= 2
        pause;
    end
    j = adjBV(1);
    k = adjBV(2);
    % detect orientation
    [ri,~] = find(F == i);
    [rj,~] = find(F(ri,:) == j);
    [rk,~] = find(F(ri,:) == k);
    Fj = F(ri(rj),:);
    Fk = F(ri(rk),:);
    ci = find(Fj == i);
    cj = find(Fj == j);
    if mod(ci+1,3) == mod(cj,3) % j is on the right hand side of i
        e1 = V(j,:) - V(i,:);
    elseif mod(cj+1,3) == mod(ci,3) % i is on the right hand side of j
        e1 = V(i,:) - V(j,:); 
    else
        pause;
    end
    ci = find(Fk == i);
    ck = find(Fk == k);
    if mod(ci+1,3) == mod(ck,3) % k is on the right hand side of i
        e2 = V(k,:) - V(i,:);
    elseif mod(ck+1,3) == mod(ci,3) % i is on the right hand side of k
        e2 = V(i,:) - V(k,:); 
    else
        pause;
    end
    Vec = -(e1 + e2);
    Vec = Vec / sqrt(sum(Vec.^2));
    Vt(idx,:) = Vec;
end
% 3. Boundary Binormal
BiN = cross(VN, Vt);
Vout = V(BV,:) ./ repmat(sqrt(sum(V(BV,:)'.^2))',[1, 3]);
iVt = BiN;
tVt = BiN;

% Construct matrices for eigs (simplified)
W = speye(4*nV,4*nV);
U = sparse(nV*4, (nV-size(BV,2))*4 + size(BV,2)*2);
U(1:(nV-size(BV,2))*4, 1:(nV-size(BV,2))*4) = speye(4*(nV-size(BV,2)), 4*(nV-size(BV,2)));
yIdx = 4*(nV-size(BV,2))+1;
tIdx = 1;
for idx = (nV-size(BV,2))+1 : nV
    U(idx*4-3:idx*4, yIdx:yIdx+1) = [1 0; 0 tVt(tIdx,1); 0 tVt(tIdx,2); 0 tVt(tIdx,3)];
    yIdx = yIdx + 2;
    tIdx = tIdx + 1;
end
U = W*U;

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

% Eigs
mat = (D*U)'*Mf*(D*U);
mat = (mat + mat') / 2;
massMat = (B*U)'*Mf*(D*U);
massMat = (massMat + massMat') / 2;
[eVectemp, eVal] = eigs(mat, massMat, 400, 'sm');
[~ ,i] = sort(diag(abs(eVal))); % sort
eVal = sum(eVal,2);
eVal = eVal(i);
eVec = U*eVectemp;
eVec = eVec(:,i);

% Dirac Disk Spin transformation
numDisk = 16;
color = ones(size(V,1),1);
for idx = [1:4:numDisk*4]
    Vnew = SpinTrans(F,V,eVec(:,idx));
    % Visualize Dirac Disk
    subplot(4,4,(idx-1)/4+1)
    trimesh(F,Vnew(:,1),Vnew(:,2),Vnew(:,3), color, 'edgealpha', 0.4,'facealpha',0.4)
    axis equal
    axis off
    shading interp 
end
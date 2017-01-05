% ================================
% Demo of Heat Kernel Signature
% This method is based on the following paper:
%
% Sun, Jian, Maks Ovsjanikov, and Leonidas Guibas. 
% "A Concise and Provably Informative Multi?Scale Signature Based on Heat Diffusion." 
% Computer graphics forum. Vol. 28. No. 5. Blackwell Publishing Ltd, 2009.
% ================================
 
clc; clear all; close all;
[V,F] = readOBJ('spot.obj');

% Construct laplace matrices
L = CotanLaplace(V, F);
[~, vertAreaList, ~, ~] = ComAreaAndNormal(F,V);
Mv = sparse(size(vertAreaList,1), size(vertAreaList,1));
for idx = 1:size(vertAreaList,1)
    Mv(idx, idx) = vertAreaList(idx);
end

% Eigs
[eVec, eVal] = eigs(L , Mv, 100, 'sm'); %+ 1e-6.*speye(size(L,1),size(L,2))
[~ ,i] = sort(diag(abs(eVal))); % sort
eVal = sum(eVal,2);

% Create time step
tmin = abs(4*log(10) / eVal(end));
tmax = abs(4*log(10) / eVal(2));
nstep = 100;
stepsize = (log(tmax) - log(tmin)) / nstep;
logts = log(tmin):stepsize:log(tmax);
ts = exp(logts);
hksOrig = eVec(:, 2:end).^2 * exp( - abs(eVal(2:end)) * ts); % original heat kernel signature

% Make it scale invariant
hks = abs(eVec(:,2:end)).^2 * exp((abs(eVal(2))-abs(eVal(2:end))) * ts);
colSum = sum(Mv*hks);
scale = 1.0./ colSum; 
scaleMat = sparse(1:length(scale), 1:length(scale), scale);
hks = hks * scaleMat; % scale invariant heat kernel signature

% Visualize HKS
for t = 90:size(hks,2)
    clf;
    trimesh(F,V(:,1),V(:,2),V(:,3),hks(:,t));
    view(120,10)
    colormap jet
    axis equal
    axis off
    pause(0.5);
end



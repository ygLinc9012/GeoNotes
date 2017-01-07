clear all; close all; clc;

% sigma - control how detailed it is for the MLS approximation 
sigma = 0.4;

%% create data points
[x, y] = meshgrid(0: 0.2: pi, 0: 0.2: pi);
z = cos(x) .* cos(2*y);
x = reshape(x, [], 1);
y = reshape(y, [], 1);
z = reshape(z, [], 1);

[xNew, yNew] = meshgrid(0.1: 0.2: pi, 0.1: 0.2: pi);
xNew = reshape(xNew, [], 1);
yNew = reshape(yNew, [], 1);
yNew = yNew - 0.1;

%% construct similarity matrix - Gaussian Linear
% compute points for compute similarity
for i = 1 : size(x,1)*size(x,2);
    xj(:,i) = x(i) * ones(size(x));
    yj(:,i) = y(i) * ones(size(y));
end
xi = x * ones(1, size(xj,2));
yi = y * ones(1, size(yj,2));

% compute distance
dist = (abs(xi-xj) + abs(yi-yj)) ./ sigma^2;
K = exp(-dist);
K = K./ (ones(size(K,1),1)*sum(K,1)); % normalize K

% calculate approximated points
zz = K * z;

%% visualize data
plot3(x, y, z, 'o');
title(['sigma = ', num2str(sigma)]);
hold on
grid on
axis equal
plot3(x, y, zz, '*');
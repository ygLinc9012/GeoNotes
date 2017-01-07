function L = CotanLaplace(V, F)
% This function Cotan-weighted Laplacian matrix for a triangular mesh
% Input:
% F: N by 3 face information
% V: N by 3 vertex positions
%
% Output:
% L: Cotan Laplacian Matrix

    numVert = size(V,1);
    numFace = size(F,1);

    temp1 = zeros(numFace,3);
    temp2 = zeros(numFace,3);
    angles = 0*F; % inner face angles

    for i=1:3
        i1 = mod(i-1,3)+1;
        i2 = mod(i  ,3)+1;
        i3 = mod(i+1,3)+1;
        temp1 = V(F(:,i2),:) - V(F(:,i1),:);
        temp2 = V(F(:,i3),:) - V(F(:,i1),:);
        % normalize the vectors
        temp1 = temp1 ./ repmat( max(sqrt(sum(temp1.^2,2)),eps), [1 3] );
        temp2 = temp2 ./ repmat( max(sqrt(sum(temp2.^2,2)),eps), [1 3] );
        % compute angles
        angles(:,i1) = acos(sum(temp1.*temp2,2));
    end

    % then compute cotan laplace L = -1/2 * (cot a + cot b) * (uj - ui) 
    L = sparse(numVert,numVert);
    for i=1:3
        i1 = mod(i-1,3)+1;
        i2 = mod(i  ,3)+1;
        i3 = mod(i+1,3)+1;
        L = L + sparse(F(:,i1),F(:,i2),-cot(angles(:,i3)),numVert,numVert,numFace);       
    end
    L = 1/2 * (L + L');
    L = sparse(1:numVert,1:numVert,-sum(L,2),numVert,numVert,numVert) + L;
end
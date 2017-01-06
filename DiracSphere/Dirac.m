function D = Dirac(F, V)
% This function eats a triangular mesh and outputs the dirac operator
%
% Inputs:
% F: N by 3 faces matrix
% V: N by 3 vertices matrix
%
% Output:
% D: quaternionic Dirac operator
%
% Reference: 
% Crane, Keenan, Ulrich Pinkall, and Peter Schröder. 
% "Spin transformations of discrete surfaces." 
% ACM Transactions on Graphics (TOG). Vol. 30. No. 4. ACM, 2011.

    [faceAreaList, ~, ~, ~] = ComAreaAndNormal(F,V);
    nV = size(V,1);
    nF = size(F,1);
    
    D = sparse(4*nF,4*nV);
    edgVec = zeros(3,4); % store 4D edge vector 
    for faceIdx = 1:nF
        face = F(faceIdx, :);
        vert = V(face, :);
        % area = norm( cross(vert(2,:)-vert(1,:), vert(3,:)-vert(1,:)) )/2;
        for idx = 1:3
            edgVec(idx,:)=[0 V(face(mod(idx+1,3)+1),:)-V(face(mod(idx+0,3)+1),:)];
        end
        for idx = 1:size(face,2)
            vertIdx = face(idx);
            edgQuatMat = quat(edgVec(idx,:));
            D(faceIdx*4-3 : faceIdx*4, vertIdx*4-3 : vertIdx*4) = ...
                D(faceIdx*4-3 : faceIdx*4, vertIdx*4-3 : vertIdx*4) + ...
                (-1/(2*faceAreaList(faceIdx))) * edgQuatMat; 
        end
    end

end
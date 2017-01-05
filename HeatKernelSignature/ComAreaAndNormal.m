function [faceArea, vertArea, faceNormal, vertNormal] = ComAreaAndNormal(F,V)
% This function compute face area, vertex area, face normal, vertex normal
% Input:
% F: N by 3 face information
% V: N by 3 vertex positions
%
% Output:
% face area, vertex area, face normal, vertex normal

    vertex = V';
    face = F';

    % Compute face normal (unnormalized)
    faceNormaltemp = cross( vertex(:,face(2,:))-vertex(:,face(1,:)), ...
        vertex(:,face(3,:))-vertex(:,face(1,:)) ) / 2;
    % normalize face normal
    faceNormal = faceNormaltemp;% ./ repmat(sqrt(sum(faceNormaltemp.^2,1)), [3 1]);
    
    % Compute face Area
    faceArea = sqrt(sum(faceNormaltemp.^2,1));
    faceArea = faceArea';
    
    vertArea = zeros(size(vertex,2), 1);
    vertNormaltemp = zeros(size(vertex,2), 3);
    for faceIdx = 1:size(face,2)
        v1 = face(1,faceIdx);
        v2 = face(2,faceIdx);
        v3 = face(3,faceIdx);
        
        % Compute vertex area
        vertArea(v1) = vertArea(v1) + faceArea(faceIdx)/3;
        vertArea(v2) = vertArea(v2) + faceArea(faceIdx)/3;
        vertArea(v3) = vertArea(v3) + faceArea(faceIdx)/3;   
        
        % Compute vertex normal(unnormalized)
        vertNormaltemp(v1,:) = vertNormaltemp(v1,:) + faceNormaltemp(:,faceIdx)';
        vertNormaltemp(v2,:) = vertNormaltemp(v2,:) + faceNormaltemp(:,faceIdx)';
        vertNormaltemp(v3,:) = vertNormaltemp(v3,:) + faceNormaltemp(:,faceIdx)';
    end
    % normalize vertex normal
    vertNormal = vertNormaltemp ./ repmat(sqrt(sum(vertNormaltemp.^2,2)), [1 3]);
end


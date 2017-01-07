function Vnew = SpinTrans(F,V,lambda)
% Build quaternionic laplace
Ltemp = CotanLaplace(V, F);
[r, c, val] = find(Ltemp);
R = [4*r-3 4*r-2 4*r-1 4*r];
R = reshape(R',[4*size(R,1),1]);
C = [4*c-3 4*c-2 4*c-1 4*c];
C = reshape(C',[4*size(C,1),1]);
val = [val val val val];
val = reshape(val',[4*size(val,1),1]);
L = sparse(R,C,val,size(V,1)*4, size(V,1)*4);

% Build quaternionic divergence matrix
F = F';
V = V';
nF = size(F,2);
nV = size(V,2);
plc = -3:0;
diver = zeros(4*nV,1);    
for i = 1:nF
    for j = 1:3
        face0 = F(mod(j-1,3)+1,i);
        face1 = F(mod(j+0,3)+1,i);
        face2 = F(mod(j+1,3)+1,i);
        vec1 = V(:,face1)-V(:,face0);
        vec2 = V(:,face2)-V(:,face0);
        cotan = dot(vec1,vec2) / norm( cross(vec1,vec2) );
        lambda1 = quat(lambda(face1*4 + plc));
        lambda2 = quat(lambda(face2*4 + plc));
        edgVec = quat([0;V(:,face2)-V(:,face1)]);
        tilda = lambda1'*edgVec*lambda1/3 + lambda1'*edgVec*lambda2/6 + lambda2'*edgVec*lambda1/6 + lambda2'*edgVec*lambda2/3;
        diver(face1*4+plc,1) = diver(face1*4+plc,1) - cotan*tilda(:,1)/2;
        diver(face2*4+plc,1) = diver(face2*4+plc,1) + cotan*tilda(:,1)/2;
    end
end   

sol = reshape(L\diver, [4 nV]);
temp = sum(sol.*sol, 1);
sol = sol / sqrt(max(temp));
Vnew = sol(2:end,:);
Vnew = Vnew';
end
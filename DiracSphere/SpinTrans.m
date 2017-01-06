function Vnew = SpinTrans(F,V,lambda)
% This function eats the triangular mesh and the Dirac eigenvector and 
% outputs the vertices of Dirac sphere
%
% Input:
% F: 3 by N face matrix
% V: 3 by M vertex matrix
% lambda: eigenvector of Dirac Operator
% 
% Output:
% Vnew: new vertex positions
% 
% This code is adapted from Dr. Keenan Crane's website found here:
% http://www.hakenberg.de/diffgeo/differential_geometry.htm#Spin

nF = size(F,2);
nV = size(V,2);
plc = -3:0;

L = sparse(4*nV,4*nV);
ome=zeros(4*nV,1);
for c1=1:nF
  for c2=1:3
    k0=F(mod(c2-1,3)+1,c1);
    k1=F(mod(c2+0,3)+1,c1);
    k2=F(mod(c2+1,3)+1,c1);
    u1=V(:,k1)-V(:,k0);
    u2=V(:,k2)-V(:,k0);
    cta=dot(u1,u2) / norm( cross(u1,u2) );
    h=jiH([cta*0.5 0 0 0]);
    ini=[k1*4+plc  k2*4+plc];
    L(ini,ini)=L(ini,ini)+[ h -h;-h h];
    if k1>k2
      k3=k1; k1=k2; k2=k3; % swap
    end
    lm1=jiH(lambda(k1*4+plc));
    lm2=jiH(lambda(k2*4+plc));
    edv=jiH([0;V(:,k2)-V(:,k1)]);
    til=lm1'*edv*lm1/3 + lm1'*edv*lm2/6 + lm2'*edv*lm1/6 + lm2'*edv*lm2/3;
    ome(k1*4+plc,1)=ome(k1*4+plc,1)-cta*til(:,1)/2;
    ome(k2*4+plc,1)=ome(k2*4+plc,1)+cta*til(:,1)/2;
  end
  if ~mod(c1,500); fprintf('.'); end
end
fprintf('\n')


ome=L\ome;
ome=reshape(ome,[4 nV]);
nrm=sum(ome.*ome,1);
ome=ome/sqrt(max(nrm));
Vnew=ome(2:end,:);

function h=jiH(pnt)
a=pnt(1); b=pnt(2); c=pnt(3); d=pnt(4);
h=[ a -b -c -d
    b  a -d  c
    c  d  a -b
    d -c  b  a];
function [Sest,d]=cosamp2(Phi,u,utrue,K)

% Cosamp algo 
% K : sparsity of Sest
% Phi : obs matrix
% u: obs vector
% Written by David Mary
% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Disclaimer: the short answer, this script is for educational purpose only.
% More Disclaimer:  http://igorcarron.googlepages.com/disclaimer

% Init.
Sest=zeros(size(utrue));
v=u;
t=1; T2=[];
while t < 101 
[k,z]=sort(abs(Phi'*v));k=flipud(k);z=flipud(z);
Omega=z(1:2*K);
T=sort(union(Omega,T2));phit=Phi(:,T);
%keyboard
b=abs(pinv(phit)*u);
[k3,z3]=sort((b));k3=flipud(k3);z3=flipud(z3);
Sest=zeros(size(utrue));
Sest(T(z3(1:K)))=abs(b(z3(1:K)));
[k2,z2]=sort(abs(Sest));k2=flipud(k2);z2=flipud(z2);
T2=z2(1:K);
v=u-Phi*Sest;
d=0;n2=norm(abs(Sest-utrue),'inf');
if n2 <1e-3
    d=1;t=1e10;
    disp('CoSaMP: success');
end
 t=t+1;
 end
 if d==0
     disp('CoSaMP: failed')
 end
 
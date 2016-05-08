function [Sest,d]=SPA(Phi,K,y,utrue)

% Subspace Pursuit Algorithm of Wai
% K : sparsity of Sest
% Phi : obs matrix
% y: obs vector
% Written by David Mary
% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Disclaimer: the short answer, this script is for educational purpose only.
% More Disclaimer:  http://igorcarron.googlepages.com/disclaimer

% Init.
[k,z]=sort(abs(Phi'*y));k=flipud(k);z=flipud(z);
T=z(1:K);phit=Phi(:,T);
yr=y-phit*pinv(phit)*y;

 t=1; 
 while t < 101    
[k2,z2]=sort(abs(Phi'*yr));k2=flipud(k2);z2=flipud(z2);  
T2=z2(1:K);
T3=sort(union(T,T2));phit3=Phi(:,T3);
xp=(pinv(phit3))*y;
[k3,z3]=sort(abs(xp));k3=flipud(k3);z3=flipud(z3);
T=T3(z3(1:K));phit4=Phi(:,T);
yr=y-phit4*(pinv(phit4))*y;
S=(pinv(phit4))*y;Sest=zeros(size(utrue));Sest(T)=abs(S);
n=norm(yr,2);
d=0;n2=norm(abs(Sest-utrue),'inf');
if n2 <1e-3
    d=1;t=1e10;
    disp('SPA: success');
end
 t=t+1;
 end
 if d==0
     disp('SPA: failed')
 end

function [d,u_new]=IRLSregcomp(phi,b,p,u_true);
% Algorithm implementing an approach suggested in
% "Iteratively Reweighted Algorithms for Compressive Sensing"
% by Rick Chartrand and Wotao Yin
% and also Chartrand & Staneva in
% "Restricted isometry properties and
% nonconvex compressed sensing"
% Written by David Mary
% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Disclaimer: the short answer, this script is for educational purpose only.
% More Disclaimer:  http://igorcarron.googlepages.com/disclaimer

epsilon = 1;
d=0; %correct reconstruction flag
u_0 = phi\b; %initialization
%u_0=(randn(size(u_0)));
u_old = u_0;
j=0; %iteration counter
m=length(b);n=length(u_true);%number of measurements/of unkowns
phi2=phi';
error=1;
while epsilon > 10^(-10)
j = j + 1;
w=((abs(u_old)).^(2) + epsilon).^(p/2-1);
v=1./w;
Q_n=diag(v,0);
u_new = Q_n * phi' * ((phi*Q_n*phi') \ b);

if lt(norm(u_new - u_old,2),epsilon^(1/2)/100)
epsilon = epsilon /10;
end
u_old = u_new;
if (j==1000)  %stop on j
    epsilon=0;
    display('IRLS : stop : j=1000 ')
end
error=norm(u_new-u_true,'inf');
if error < 1e-3
    epsilon=0;
    disp('IRLS : success ')
    d=1; break
end
end
% d : correct reconstuction flag
if d==0;
    disp('IRLS : failed')
end   

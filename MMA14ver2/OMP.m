function [Sest,d] =OMP(Phi,y,strue,K,tol1);

% Othogonal matching pursuit of Tropp et Gilbert
% Find greedy approx of s s.t. phi*s = y
% y : data
% phi : measurement matrix
% K :number of vectors used in the approx (typ. = S or M)
% Written by David Mary
% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Disclaimer: the short answer, this script is for educational purpose only.
% More Disclaimer:  http://igorcarron.googlepages.com/disclaimer

Sest=zeros(size(strue));
Pos=[]; % positions indexes of components of s
rtm1=y; % first residual
Phit=[]; % Matrix of the columns used to represent y 
t=1; %iteration number
d=0; % error flag
%for t=1:length(Phi(:,1));
    for t=1:K;
    [i,j]=max(abs(Phi'*rtm1));
    Pos=[Pos j];
    Phit=[Phit Phi(:,j)];
    sest=Phit\y;
    rtm1=y-Phit*sest;
    Sest(Pos)=sest;
error=norm(abs(Sest-strue),'inf');
if error < tol1
    d=1;
    disp('OMP : success')
    break
end
end;
if d==0
    disp('OMP : failed')
end

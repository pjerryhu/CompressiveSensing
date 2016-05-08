function [Sest,d] =OMPmod3(Phi,y,strue,K);

% Othogonal matching pursuit of Tropp et Gilbert
% modified  : finds at each iteration the vector
% that maximizes the loss in the residual
% and also makes sure the selected vectors have
% no linear dependencies
% Implementation probably computationnally very inefficient
% Finds greedy approx of s s.t. phi*s = y
% y : data
% phi : measurement matrix
% K :number of vectors used in the approx (typ. = S or M)
% Written by David Mary
% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Disclaimer: the short answer, this script is for educational purpose only.
% More Disclaimer:  http://igorcarron.googlepages.com/disclaimer
%
Sest=zeros(size(strue));
Pos=[]; % positions indexes of components of s
rtm1=y; % first residual
Phit=[]; % Matrix of the colums used to represent y 
t=1; %iteration number
d=0; % error flag
Mat=Phi;
    for t=1:K;
        for v=1:length(Mat(1,:))
            Phit2=[Phit Mat(:,v)]; % takes every new vector
            
            if rank(Phit2)~=length(Phit2(1,:))
                score(v)=1e20;
            else
                sest=Phit2\y;r=y-Phit2*sest;
            score(v)=norm(r,2);
            end
        end     
        [i,j]=min(score);
         Phit=[Phit Mat(:,j)]; % add selected vector to Phit
         % remove now the selected vector to Mat
       if j==1
    Mat=Mat(:,2:length(Mat(1,:)));
end
if j==length(Mat(1,:))
    Mat=Mat(:,1:length(Mat(1,:))-1);
end
if (j~=1) & (j~=(-1))
    Mat=[Mat(:,1:j-1) Mat(:,j+1:length(Mat(1,:)))];
end
    sest=Phit\y;
    rtm1=y-Phit*sest;
    % finds what is the position in Phi of the jth vector in Mat
    if t==1
     Pos=[Pos j];
    else
        l=length(find(Pos<=j));t=j+l;
        while length(find(Pos<=t))~=l
            t=t+1;l=l+1;
        end
        Pos=[Pos t];
    end
    Sest(Pos)=sest;
error=norm(abs(Sest-strue),'inf');
if error < 1e-3
    d=1;
    disp('modified OMP : success')
    break
end
%t=length(Phi(1,:));
end;
if d==0
    disp('OMP mod : failed')
end
  
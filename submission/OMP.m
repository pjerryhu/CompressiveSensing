function [Sest,d] =OMP(Phi,y,strue,K);

% Othogonal matching pursuit of Tropp et Gilbert
% Find greedy approx of s s.t. phi*s = y
% y : data
% phi : measurement matrix
% K :number of vectors used in the approx (typ. = S or M)
%    also the number of iteration

Sest=zeros(size(strue));
Pos=[]; % positions indexes of components of s
rtm1=y; % first residual
Phit=[]; % Matrix of the columns used to represent y 
t=1; %iteration number
d=0; % error flag
%for t=1:length(Phi(:,1));


for t=1:K
    [i,j]=max(abs(Phi'*rtm1));  % find the index of max term 
                                % that the easy optimization problem
    Pos=[Pos j];                % include the column into Pos matrix
    Phit=[Phit Phi(:,j)];       
    sest=Phit\y;
    rtm1=y-Phit*sest;           % subtract the largest component from estimation
    Sest(Pos)=sest;             % update the size of residual
    error=norm(abs(Sest-strue),'inf'); 
    
%     figure;
%     plot(1:length(strue),strue, 'b*')
%     hold on
%     plot(1:length(Sest),Sest,'ro');
    
%     pause(0.5)
    if error < 1e-3
        d=1;
%         disp('OMP : success')
        break
    end
end;
if d==0
%     disp('OMP : failed')
end

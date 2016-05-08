
% - OMP: this script is an implementation of the classical Orthogonal 
% Matching Pursuit algorithm.


close all
clear
tic
l=8;
c=8; 

% (l,c) : dimension of the image
% m is the number of Compressed Sensing measurements
% p is the defining the order of the metric p=1 is l1 norm....
% S is the support of the sparse signal
S=5;
% m is the number of measurement
%m=20;
A=5:20;
D=zeros(size(A));D2=D;%flags of correct reconstructions for Monte Carlo
% for t=1:5;% loop total number of experiments
% for k=1:length(D);% loop of varying figure : either S or M
k = 1
m=A(k); %assumes S was fixed above
  
 % Unknown image u_true   
 r2=randperm(l*c);r2=r2(1:S);rr2=r2(:);
 u_true=zeros(l*c,1); u_true(r2)=u_true(r2)+rand(S,1);  
 

%CS matrix phi :
% phi=(kron(dftmtx(c).',dftmtx(l)));
phi = (kron((fft(eye(c))).',fft(eye(l))));
r = randperm(l*c);r=r(1:m);
phi=phi(r,:); % if subsampled Fourier 2D matrix
%phi=randn(m,l*c); % if Gaussian matrix

b = phi * u_true(:); % observation vector
%b_noise=.001*randn(size(u_true))+j*0.001*randn(size(u_true)); %if noise
%b = b+b_noise; %vectorized image

           
[uestOMP,d2] =OMP(phi,b,u_true,S);
           
D2(k)=D2(k)+d2;


% OMP: this script is an implementation of the classical Orthogonal 
% Matching Pursuit algorithm.


close all
clc
clear

% (l,c) : dimension of the image
% m is the number of Compressed Sensing measurements
% p is the defining the order of the metric p=1 is l1 norm....
% S is the support of the sparse signal
S=5;
l = 8;
c = 8;
m = 5
  
 % Unknown image u_true   
 r2=randperm(64);r2=r2(1:5);rr2=r2(:);
 u_true=zeros(64,1); u_true(r2)=u_true(r2)+rand(S,1);  
 

%CS matrix phi :
phi = randn(15,64);
% phi = (kron((fft(eye(c))).',fft(eye(l))));
% r = randperm(l*c);r=r(1:m);
% phi=phi(r,:); % if subsampled Fourier 2D matrix


b = phi * u_true; % observation vector
%b_noise=.001*randn(size(u_true))+j*0.001*randn(size(u_true)); %if noise
%b = b+b_noise; %vectorized image

           
[uestOMP,d2] =OMP(phi,b,u_true,9);



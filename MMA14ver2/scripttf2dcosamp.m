% This script checks the performance of some reconstruction
% used in Compressed Sensing / Compressive Sensing. It includes the
% following implementations:
%
% - IRLSregcomp : This script implements an adaptation
% of the algorithm of Rick Chartrand and Valentina Staneva featured in 
% "Restricted isometry properties and nonconvex compressed sensing"
%
% - OMP: this script is an implementation of the classical Orthogonal 
% Matching Pursuit algorithm.
%
% - OMPmod3 : this script is a modified version of OMP. It is a little more
% but it is also more efficient
%
% - CoSaMP : this script is an implementation of CoSaMP from
%"CoSaMP: Iterative signal recovery from incomplete and inaccurate samples"
% by Deanna Needell and Joel Tropp.
%
% - SPA : An implementation of "Subspace Pursuit for Compressive Sensing:
% Closing the Gap Between Performance and Complexity" by Wei Dai and 
% Olgica Milenkovic.
%
% The original image to be reconstructed is made of pixels 
% of various intensities laying on a dark background (think of stars 
% at night). It is sparse in the direct/image space.
%
% The measurement matrix is a 2D Random Partial Fourier Transform.
%
% The goal is t find a sparse approximation of u_true
% knowing that b=phi*u_true (+ noise)
% b: obs vector, phi : CS matrix, u_true : target object
%
% Written by David Mary
% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Disclaimer:the short answer, this script is for educational purpose only.
% More Disclaimer:  http://igorcarron.googlepages.com/disclaimer
%

close all
clear
tic
l=8;
c=8; 
tol1 = 1e-8;

% (l,c) : dimension of the image
% m is the number of Compressed Sensing measurements
% p is the defining the order of the metric p=1 is l1 norm....
% S is the support of the sparse signal
S=5;
% m is the number of measurement
%m=20;
A=5:20;
D=zeros(size(A));D2=D;Dmod=D;D1=D;Dp5=D;Dmod3=D;Dco=D;Dsp=D; %flags of correct reconstructions for Monte Carlo
for t=1:200;% loop total number of experiments
for k=1:length(D);% loop of varying figure : either S or M
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

           [d,u_new]=IRLSregcomp(phi,b,0.01,u_true(:),tol1);
           [dp5,u_newp5]=IRLSregcomp(phi,b,0.5,u_true(:),tol1);
           [d1,u_new1]=IRLSregcomp(phi,b,1,u_true(:),tol1);
           [uestOMP,d2] =OMP(phi,b,u_true,S,tol1);
           [uestOMPmod3,dmod3] =OMPmod3(phi,b,u_true,S,tol1);
           [uestco,dco] =cosamp2(phi,b,u_true,S,tol1);
           [uestsp,dsp] =SPA(phi,S,b,u_true,tol1);
           
           Dsp(k)=Dsp(k)+dsp;
           Dco(k)=Dco(k)+dco;
           D(k)=D(k)+d;
           Dmod3(k)=Dmod3(k)+dmod3;
           D2(k)=D2(k)+d2;
           D1(k)=D1(k)+d1;
           Dp5(k)=Dp5(k)+dp5;
           
end
end;
           D=D/t;Dp5=Dp5/t;D1=D1/t;Dco=Dco/t;Dsp=Dsp/t;
           Dmod=Dmod/t;
           Dmod3=Dmod3/t;
           D2=D2/t;
           figure(1)
           hold off
           plot(A,D,'g*--','LineWidth',2)
           hold on
           plot(A,Dp5,'c*--','LineWidth',2)
           plot(A,D1,'k*--','LineWidth',2)
           plot(A,D2,'m*--','LineWidth',2)
           plot(A,Dmod3,'r*--','LineWidth',2)
            plot(A,Dsp,'*--','LineWidth',2)
             plot(A,Dco,'o--','LineWidth',2)
           legend('IRLSL p=.01','IRLS p=.5','IRLS p=1','OMP','modified OMP','SP','CoSaMP',4)
           grid
           xlabel('Number of measurements M')
           ylabel(['Empirical probability of reconstruction at level ' num2str(tol1)])
           %title('N=64 S=5 \Phi : undersampled DFT')
         
% u_new=reshape(u_new,l,c);
% u_true=reshape(u_true,l,c);
% figure(1)
% hold off
% imshow(real(u_true),'InitialMagnification','fit'); colormap(jet); colorbar
% title('original')
% figure(2)
% hold off
% imshow(real(u_new),'InitialMagnification','fit'); colormap(jet); colorbar
% title('reconstructed')
% figure(3)
% hold off
% imshow(u_true-real(u_new),'InitialMagnification','fit'); colormap(jet); colorbar
% title('error')
% toc
% % locations of the points in the TF
% figure(4)
% hold off
% loc=zeros(l*c,1);loc(r)=ones(size(r));loc=reshape(loc,l,c);
% imshow(loc,'InitialMagnification','fit')
% title('(u,v) coverage')
% OMP: this script is an implementation of the classical Orthogonal 
% Matching Pursuit algorithm.

% oldErrors = errors
close all
clc
clear

% (l,c) : dimension of the image
% m is the number of Compressed Sensing measurements
% p is the defining the order of the metric p=1 is l1 norm....
% S is the support of the sparse signal
figure1 = figure
hold on
start = toc;
for S=5:5:40

cols = 64;
tic;
% Unknown image u_true   
r2=randperm(cols);r2=r2(1:S);rr2=r2(:);
u_true=zeros(cols,1); u_true(r2)=u_true(r2)+randn(S,1);  


%CS matrix phi :

% phi = (kron((fft(eye(c))).',fft(eye(l))));
% r = randperm(l*c);r=r(1:m);
% phi=phi(r,:); % if subsampled Fourier 2D matrix

K = 4;
t = 50;
% for test = 1:5
% for sigma = 0.01:0.01:0.36
sigma = 0.36;
Nmax = ceil(K*S*log(cols/sigma));
if Nmax >= cols
    Nmax = cols;
end
%%
errors = zeros(t,length(5:5:ceil(Nmax/5)*5));
stride = 5;
for test = 1:t
for N = stride:stride:ceil(Nmax/stride)*stride
    
    A = randn(N,cols);
    b = A * u_true; % observation vector
    %b_noise=.001*randn(size(u_true))+j*0.001*randn(size(u_true)); %if noise
    %b = b+b_noise; %vectorized image
    [uestOMP,d2] = OMP(A,b,u_true,N);
    
    errors(test,N/stride) = norm(uestOMP-u_true,2);
    
end
end
endd = toc

%%

errormean = mean(errors);

N = stride:stride:ceil(Nmax/stride)*stride;
if S == 5
    plot(N,errormean,'b-o')
elseif  S == 10
    plot(N,errormean,'r-o')
elseif  S == 15
    plot(N,errormean,'g-o')
elseif  S == 20
    plot(N,errormean,'b-*')
elseif  S == 25
    plot(N,errormean,'r-*')
elseif  S == 30
    plot(N,errormean,'g-*')
elseif  S == 35
    plot(N,errormean,'b-^')
elseif  S == 40
    plot(N,errormean,'r-^')
end
end
%%
legend('S = 5','S = 10','S = 15','S = 20','S = 25','S = 30','S = 35','S = 40')
title('Monte Carlo Simulation of OMP error rate with increasing measurements from 5 to c*S*log(cols/delta)')
xlabel('N: # of measurements');ylabel('deviation: ell-2 norm of recovered signal - actual signal');
saveas(figure1,'OMP_withNrun.png')

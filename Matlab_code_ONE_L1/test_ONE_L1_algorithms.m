%%
%%%%%%%%%%%%%%%%%%% rONE-L1 %%%%%%%%%%%%%%%%%%%%
% random DCT matrix ensemble and real-valued signal
clear
clc
close all
N = 2^12;     % data dimension
delta = .3;   % sampling ratio  n/N
rho = .2;    % sparsity      k/n

n = ceil(delta*N);
k = ceil(n*rho);

p1 = randperm(N);
pos = p1(1:k);
x0 = zeros(N,1);
x0(pos) = randn(k,1);

p2 = randperm(N);
Omega = p2(1:n)';

A = pDCT(N,Omega);

b = A*x0;

% r = 1.1;
% r = 1.02;
r = min(1+0.04*delta, 1.02);
tol = 1e-5;
maxiter = 1000;

tstart = tic;
% [x,iter] = greedyIST(A,b,r,tol,maxiter,0);
[x, iter] = rONE_L1(A, b,r,tol,maxiter,0);
time = toc(tstart);
dif = norm(x-x0,2)/norm(x0,2);

fprintf('=== Reconstruction using rONE-L1 ===\n');
fprintf('Real-valued settings: random DCT matrix emsemble and real signal\n\n');
fprintf('Relative root mean squared error: %e \nNo. of iterations: %d\nCPU time: %f\n\n\n\n', dif, iter, time);


%%
%%%%%%%%%%%%%%%%%% explicit matrix %%%%%%%%%%%%%%%%%%%
clear all

N = 1000;     % data dimension
delta = .3;   % sample ratio  n/N
rho = .2;    % sparsity      k/n

n = ceil(delta*N);
k = ceil(n*rho);

p1 = randperm(N);
pos = p1(1:k);
x0 = zeros(N,1);
x0(pos) = randn(k,1);

A1 = randn(n,N);

A = orth(A1.').';
b = A*x0;
r = min(1+0.04*delta, 1.02);
tol = 1e-5;
maxiter = 1000;

tstart = tic;
[x, iter] = rONE_L1(A,b,r,tol,maxiter,0);
time = toc(tstart);

dif = norm(x-x0,2)/norm(x0,2);

fprintf('=== Reconstruction using rONE-L1 ===\n');
fprintf('Real-valued settings: explicit matrix and real signal\n\n');
fprintf('Relative root mean squared error: %e \nNo. of iterations: %d\nCPU time: %f\n\n\n\n', dif, iter, time);


%%
%%%%%%%%%%%% random FFT matrix ensemble, complex-valued signal %%%%%%%%%%%
clear all

N = 2^12;     % data dimension
delta = .3;   % sample ratio  n/N
rho = .2;    % sparsity      k/n

n = ceil(delta*N);
k = ceil(n*rho);

p1 = randperm(N);
pos = p1(1:k);
x0 = zeros(N,1);
% x0(pos) = randn(k,1);
x0(pos) = randn(k,1) + 1i*randn(k,1);

p2 = randperm(N);
Omega = p2(1:n)';

A = pDFT(N,Omega);

b = A*x0;

% r = 1+delta;
% r = 1.1;
r = min(1+0.04*delta, 1.02);
tol = 1e-5;
maxiter = 1000;

tstart = tic;
[x, iter] = rONE_L1(A,b,r,tol,maxiter,0);
time = toc(tstart);

dif = norm(x-x0,2)/norm(x0,2);

fprintf('=== Reconstruction using rONE-L1 ===\n');
fprintf('Complex-valued settings: random FFT matrix ensemble and complex signal\n\n');
fprintf('Relative root mean squared error: %e \nNo. of iterations: %d\nCPU time: %f\n\n\n\n', dif, iter, time);




%%
%%%%%%%%%%%%%%% eONE-L1 %%%%%%%%%%%%%%%%%%%%%%%
% random FFT matrix ensemble, complex-valued signal
clear all

N = 2^12;     % data dimension
delta = .3;   % sample ratio  n/N
rho = .2;    % sparsity      k/n

n = ceil(delta*N);
k = ceil(n*rho);

p1 = randperm(N);
pos = p1(1:k);
x0 = zeros(N,1);
% x0(pos) = randn(k,1);
x0(pos) = randn(k,1) + 1i*randn(k,1);

p2 = randperm(N);
Omega = p2(1:n)';

A = pDFT(N,Omega);

b = A*x0;

r = 1+delta;

tol = [1e-5 1e-6];
maxiter = [500 500];

tstart = tic;
[x, iter_tot, iter] = eONE_L1(A,b,r,tol,maxiter,0);
time = toc(tstart);

dif = norm(x-x0,2)/norm(x0,2);

fprintf('=== Reconstruction using eONE-L1 ===\n');
fprintf('Complex-valued settings: random FFT matrix ensemble and complex signal\n\n');
fprintf('Relative root mean squared error: %e \nTotal No. of iterations(inner loop): %d\n', dif, iter_tot);
fprintf('No. of iterations(outer loop): %d\nCPU time: %f\n\n', iter, time);


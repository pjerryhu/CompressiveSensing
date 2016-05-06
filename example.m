%%
close all
clc
clear

n = 512;
m = 256;

A = randn(m,n);
S = round(m/5);

support = randsample(n,S);
x0 = zeros(n,1); x0(support) = randn(S,1);
b = A*x0;

%%

cvx_begin

variable x(n);
minimize(norm(x0,1));

subject to
    A*x == b;

cvx_end

norm(x-x0)/norm(x0)
figure; 
plot(1:n, x0, 'b*', 1:n, x, 'ro')
%% Basis Pursuit w/o white noise
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

plot(1:n, x0, 'b*', 1:n, x, 'ro', 1:n, result,'go')



%% Basis Pursuit with white noise
close all
clc
clear

p = 512;
n = 256;

A = randn(n,p)/sqrt(n);
S = round(n/4);

support = randsample(p,S);
x0 = zeros(p,1); x0(support) = randn(S,1);
Ax0 = A*x0;
epislon = norm(Ax0)/S; % this is the signal to noise ratio
z = randn(n,1);

y = Ax0 + z*epislon/norm(z);

%%

cvx_begin

variable x(p);
minimize(norm(y,1));

subject to
     norm(A*x - y,2)<= epislon;

cvx_end

norm(x-x0)/epislon
figure; 
plot(1:p, x0, 'b*', 1:p, x, 'ro')
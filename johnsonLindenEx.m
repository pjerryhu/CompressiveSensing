clear
% Fast Johnson-Lindenstrauss Transform
% as described by Nir Ailon and Bernard Chazelle
% in their paper:
% "Approximate Nearest Neighbors and the Fast Johnson-Lindenstrauss
% Transform"
% The algorithm uses the PHD transform of paragraph 2
% The transform is a k x d projection
% for every set S of n points in R^d
% Let us define an example case
% All columns of G are sets of points
% i.e. we have 1000 points of dimension 128
n = 1000;
d = 256;
% let us examine a reduction from d = 256 to k = 150
%
k = 150;
G = 180 * rand(d,n);
% k = O(epsilon^(-2) * log(n))
epsilon = sqrt(log(n)/k)
%
% epsilon is really O( sqrt(log(n)/k) )
%
%
% Projection in dim k << d % let us take d = 2^w w = log2(d); % % % Defining P (k x d) % first define q % q = min((log(n))^2/d,1); % P = rand(k,d); P = (P < q); P1 = 1/q * randn(k,d); P = P1 .* P; % % Defining H (d x d) H = [1]; for i=1:w H = [ H H;H -H]; end H = 1/d^(0.5)* H; % % Defining D ( d x d) vec1=rand(d,1); v1 = ((vec1 > 0.5) - 1*(vec1<=0.5)); D = diag(v1); % % % FJLT = P * H * D; u = FJLT * G(:,5); v = FJLT * G(:,36); norm(G(:,5)-G(:,36))^2 * k * (1-epsilon) norm( u - v)^2 norm(G(:,5)-G(:,36))^2 * k * (1+epsilon)

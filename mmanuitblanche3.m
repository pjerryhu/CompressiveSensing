clear all;
clc
close all
% Algorithm implementing Recognition via Sparse Representation
% or Algorithm 1 suggested in
% "Feature selection in face recognition: A sparse representation
% perspective" by Allen Y. Yang, John Wright, Yi Ma, and Shankar Sastry
% written by Jort Florent Gemmeke and Igor Carron.
  
% n total number of images in the training database
n = 20;
% m is dimension of the image/ambient 
% space (number of samples in signal ) 
% (e.g. total "training" set) ..e.g. 101
x=[0:0.01:1];
m = size(x,2);  
% d is the number of Compressed Sensing 
% measurements, it is also the dimension of the feature 
% space. 
d = 7;     
% k is the number of classification groups and  for k subjects
k = 10;   
% e.g. every group consists of kn = 20 samples.
kn = n/k; 

% Building the training database/dictionary A of n elements
% Within that dictionary there are k classes
% Each of these classes correspond to functions:cos(i x)
% with added noise. i is the differentiator between classes.
% Within each class, any element is different from the other
% by some amount of noise.
% Each element is listed in a row. Each function-element is 
% sampled m times.
noise=0.1;
for i=1:k
   for j=1:kn
      yy = noise*rand(1,m);
      A(:,(i-1)*kn+j)=cos(i.*x)+yy;
   end
end   

% In the following, every index i contains a vector
% of length m with nonzeros corresponding to 
% columns in A belonging to class i (i in [1,k])
% Every row holds one element of a class, 
% i.e. there are as many rows as elements in this 
% dictionary/database.
for i=1:k
   onesvecs(i,:)=[zeros(1,(i-1)*kn) ones(1,kn) zeros(1,n-((i-1)* kn + kn))];
end
% This is the element for which we want to 
% know to what class it belongs to. Here we give
% the answer for the plot in x_true.
x_true = [zeros(1,4) 1 1 zeros(1,n-6)];   
x_true = x_true/norm(x_true,1);
% Now the expression stating that y = cos(3.*x);
y = cos(3.*x)';
% Please note it is different from both entries
% in the dictionary as there is no noise.
% we could have put this as well: y = A * x_true'
% but it would have meant we were using the training
% function as a query test.
%
%
%       
% Step 2
% Dimensionality reduction from 101 to 7 features 
% using random projection.
R = randn(d,m);         
Atilde = R * A;
ytilde = R * y;

% Normalize columns of Atilde and ytilde
for i=1:size(A,2)
   Atilde(:,i)=Atilde(:,i)/norm(Atilde(:,i),2);
end
ytilde = ytilde/norm(ytilde);

% Step 3
% call any L1 solver you want. Here we call GPSR 
%
% ---- L1 Magic Call
%xp = l1qc_logbarrier(Atilde\ytilde, Atilde, [], ytilde, 0.01, 1e-3);
%
% ---- Sparsify Solver Call
%xp = greed_omp(ytilde,Atilde,n);      
%
% ---- CVX Solver Call
%sigma = 0.05;
%cvx_begin
%	variable xp(k); 
%	minimize(norm(xp,1));
%	subject to
%		norm(Atilde*xp - ytilde) < sigma;
%cvx_end
%
% ---- GPSR Solver call
xp = GPSR_BB(ytilde,Atilde, 0.01*max(abs(Atilde'*ytilde)),'Continuation',1);

% plot results of L1 solution
figure; plot(xp,'v');
title('L1 solution (v), initial signal (o)');
hold; plot(x_true,'o');

%
% Step 4
% Calculate residuals (formula 13)
residual = [];
for i=1:k
   % only keeps nonzeros in xp where onesvec is nonzero
   deltavec(:,i) = onesvecs(i,:)'.* xp;         
   % calculate residual on only a subset of Atilde
   residual(i) = norm(ytilde-Atilde*deltavec(:,i));      
end
figure; bar(residual);                     % plot residuals
title('The lowest bar indicates the class number to which the element belongs to')

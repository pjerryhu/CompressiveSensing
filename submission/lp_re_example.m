
close all
clc
clear

figure1 = figure
hold on
p = 1;
cols = 64;
t = 50;
Nmax = cols;

tic
for S=5:5:40

% Unknown image u_true   
r2=randperm(cols);r2=r2(1:S);rr2=r2(:);
u_true=zeros(cols,1); u_true(r2)=u_true(r2)+randn(S,1);  




%%
errors = zeros(t,length(5:5:floor(Nmax/5)*5));
stride = 5;
for test = 1:t
for N = stride:stride:floor(Nmax/stride)*stride
    
    A = randn(N,cols);
    b = A * u_true; % observation vector
    uestOMP = lp_re(A,b,p);
    
    errors(test,N/stride) = norm(uestOMP-u_true,2);
    
end
end


%%

errormean = mean(errors);

N = stride:stride:floor(Nmax/stride)*stride;
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
endd = toc
legend('S = 5','S = 10','S = 15','S = 20','S = 25','S = 30','S = 35','S = 40')
title('Monte Carlo Simulation of Non-convex error rate with increasing measurements from 5 to 60')
xlabel('m: # of measurements');ylabel('deviation: ell-2 norm of recovered signal - actual signal');
saveas(figure1,'LP_RE_withNrunP1.png')

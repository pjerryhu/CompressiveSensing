function x = BasisPursuit(A,u_true)

b = A * u_true;
[m,n]=size(A);
cvx_begin

variable x(n);
minimize(norm(u_true,1));

subject to
    A*x == b;

cvx_end



end
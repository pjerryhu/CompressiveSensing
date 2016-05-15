function x = lp_re(A,y,p)
%
% lp-Reweighted Least Squares Problem Solver
% Algorithm implementing an approach suggested in
% "Iteratively Reweighted Algorithms for Compressive Sensing"
% by Rick Chartrand and Wotao Yin
% Algorithm implemented as featured in:
% http://math.lanl.gov/Research/Publications/Docs/chartrand-2008-iteratively.pdf
%
%
%   lp_re solves problems of the following form:
%				minimize ||x||_p
%
%           subject to A*x = y 
%
%   where A and y are problem data and x is variable (described below).
%
% CALLING SEQUENCES
%   [x] = lp_re(A,y,p);
%
% INPUT
%   A       : mxn matrix; input data. columns correspond to features.
%
%   m       : number of measurements ( or rows) of A
%   n       : dimension of ambient space (or columns) of A
%
%   y       : vector of result.
%
%   p       : norm order; p = 1 is l1 norm
%
% OUTPUT
%   x       : vector of vector n; initial unknown
%
epsilon = 1;
% u_0 is the L_2 solution which would be exact if m = n,
% but in Compressed expactations are that m is less than n
u_0 = A\y;
u_old = u_0;
j=0;
while epsilon > 0.00001
	j = j + 1;
	w = (u_old.^(2) + epsilon).^(p/2-1);
	v = 1./w;
	Q_n = diag(v,0);
	tu = inv(A*Q_n*A');
	u_new = Q_n * A' * tu * y;
	if lt(norm(u_new - u_old,2),epsilon^(1/2)/100)
		epsilon = epsilon /10;
	end
	u_old = u_new;
end
x=u_new;
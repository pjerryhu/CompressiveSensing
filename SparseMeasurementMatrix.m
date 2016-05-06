function [A]=SparseMeasurementMatrix(m,n,d)
% Algorithm implementing the Sparse Measurement Matrix
% construction of Radu Berinde and Piotr Indyk in 
% "Sparse Recovery Using Sparse Random Matrices" 
% MIT-CSAIL-TR-2008-001
%  available at: http://hdl.handle.net/1721.1/40089
%
% INPUT
%
% m 			: 	number of rows (Number of measurements)
%              generally m should be much less than n
%              n is limited by (m,d) or 
%              factorial(m)/(factorial(d)*factorial(m-d))
%
% n 			: 	number of columns (signal ambient dimension)
%
% d 			: 	number of ones per column
%
%
% CALLING SEQUENCES
%
%   A=SparseMeasurementMatrix(m,n,d);
%
% OUTPUT
%   A       :  sparse mesurement matrix.
%
% USAGE EXAMPLES
%   
%	A=SparseMeasurementMatrix(m,n,d);
%
% AUTHOR    Igor Carron <igorcarron@gmail.com>
% UPDATE    ver. 0.1, Feb 7 2008
%
% Creative Commons Licence
% http://en.wikipedia.org/wiki/Creative_Commons
%
A=zeros(m,n);
ntotal = nchoosek(m,d); 
if n > ntotal
   sprintf('%s','Not enough measurements, sorry')
   return
end   
jj = 0;
for i=1:n
	while sum(A(:,i),1)~=d | jj==47
		jj=0;
		A(:,i)=zeros(m,1);
		for j=1:d
			A(round((m-1)*rand(1,1))+1,i)=1;
		end
		v1 = diag(A(:,i))'*ones(m,n);
		w = A - v1;
		for j=1:i-1
			if w(:,j)== zeros(m,1)
				jj = 47;
			end
		end
	end
end
% Measurement Matrix is A
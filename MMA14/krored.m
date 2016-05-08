function [phi]=krored(A,B,r);

% A and B : squared matrices of size (NxN)
% r : vecteur m x 1 of integers correspondig to the rows
% of kron(A,B) to be selected 
% The matrix with the selected rows is phi and
% corresponds to a mask obtained in interferometry
% Written by David Mary
% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Disclaimer: the short answer, this script is for educational purpose only.
% More Disclaimer:  http://igorcarron.googlepages.com/disclaimer

N=length(A(:,1));

for j=1:length(r);
    k=floor((r(j)-1)/N);
    s=r(j)-k*N;
    v=[];
    for p=1:N;
        v=[v A(k+1,p)*B(s,:)];
    end
    phi(j,:)=v;
end
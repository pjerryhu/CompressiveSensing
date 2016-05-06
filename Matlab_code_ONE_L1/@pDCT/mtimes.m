function res = mtimes(a,b)

if a.adjoint
    bfull = zeros(a.len,1);
    bfull(a.Omega) = b;
%     res = ifft(bfull)*sqrt(a.len);
    res = idct(bfull);
else
%     fb = fft(b)/sqrt(a.len);
    fb = dct(b);
    res = fb(a.Omega);
end



    

function [ x ] = trin( L,b )
% Paola Mateos 118798
% Gabriel S�nchez 120141
% Alejandro �lvarez 113638
% Soluci�n del sistema lineal; L*x = b
% In
% L: triangular inferior no singular nxn.
% b: vector de dimensi�n n.
% Out
% x: un vector de dimensi�n n, soluci�n al sistema lineal.

n = length(b);
x = zeros(n,1);
for k = 1:(n-1)
    x(k) = b(k)/L(k,k);
    b(k+1:n) = b(k+1:n)-x(k)*L(k+1:n,k);
end
x(n) = b(n)/L(n,n);
end


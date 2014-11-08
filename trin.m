function [ x ] = trin( L,b )
% Paola Mateos 118798
% Gabriel Sánchez 120141
% Alejandro Álvarez 113638
% Solución del sistema lineal; L*x = b
% In
% L: triangular inferior no singular nxn.
% b: vector de dimensión n.
% Out
% x: un vector de dimensión n, solución al sistema lineal.

n = length(b);
x = zeros(n,1);
for k = 1:(n-1)
    x(k) = b(k)/L(k,k);
    b(k+1:n) = b(k+1:n)-x(k)*L(k+1:n,k);
end
x(n) = b(n)/L(n,n);
end


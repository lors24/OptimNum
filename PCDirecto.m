function [ x, lambda ] = PCDirecto(Q, A, c, b)
% 22/08/14
% Optimizaci?n numerica
% 
% Resolvemos el problema cuadr?tico
%   Min (1/2)*x'*Q*x + c'*x
%    s.a.   A*x = b
% por factorizacion LU

% In
% Q .- matriz sim?trica positiva definida de orden n
% A .- matriz mxn con rango(A) = m ( m <= n)
% c.- vector columna de dimensi{on n.
% b.- vector columna de dimensi?n m.

n = length(c);
m = length(b);

K = [ Q A'; A zeros(m)];
ld = [- c ; b];

% Factorizacion LU de K
[L, U] = lu(K);
% Sistema triangular inferior L*y = ld
y = L\ld;
% Sistema triangular superiro U*w = y
w = U\y;

x = w(1:n);
lambda = y(n-1:n+m);
end


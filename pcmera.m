function [x, y ] = pcmera(Q, A, c, b)
% Programaci?n cuadr?tica del m?todo del rango
%
% Resolvemos el problema cuadr?tico
%   Min (1/2)*x'*Q*x + c'*x
%    s.a.   A*x = b
% por el m?todo del rango

% In
% Q .- matriz sim?trica positiva definida de orden n
% A .- matriz mxn con rango(A) = m ( m <= n)
% c.- vector columna de dimensi{on n.
% b.- vector columna de dimensi?n m.

% Out

Qinv = inv(Q);
B = A*Qinv*A';
C = -A*Qinv*c-b;

R = chol(B);
y = R'\C;
lambda = R\y;

ld2 = -c - A'*lambda;
R2 = chol(Q);
y1 = R2'\ld2;
x = R2\y1;

end


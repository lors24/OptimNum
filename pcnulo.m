function [ x] = pcnulo( Q, A, c, b)
% M?todo del espacio nulo para
% Min (1/2)*x'*q*x + c'*x
% S.A. A*x = b
% donde Q es una matrix nxn, A es mxn con m<= n
% rango(A) = m y Q es sim?trica positiva definida.

Z = null(A); % base ortonormal del espacio nulo de A.
xp = A\b;    % soluci?n particular de A*x = b
Q1 = Z'*Q*Z; % matriz sim?trica positiva definida en R^(n-m)
c1 = Z'*(Q*xp + c); % lado derecho
% resolver el sistema lineal: Q1*xz = -c1
% Cholesky

% despues vamos a sustituir esto por gradiente conjugado

L = chol(Q1)';
%y = trin(L, -c1);   %Mi rutina para sistemas triangulares inferiores
%xz = tris(L',y);   %Mi rutina para sistemas triangulares superiores

y = -L\c1;
xz = L'\y;

%
x = xp + Z*xz; %solucion final




end


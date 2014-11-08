function   x = tris(U,b)
% Paola Mateos 118798
% Gabriel Sánchez 120141
% Alejandro Álvarez 113638
% Se resuelve el sistema triangular inferior Ux = b usando sustitución
% hacia atrás
% In
% U: Matriz triangular superior.
% b: vector de tamaño n, correspondiente al lado derecho de la ecuación del
% sistema.
%
% Out
% x: vector solución de tamaño n.

n = length(b);
x = zeros(n,1);

x(n) = b(n)/U(n,n);
for k = n-1:-1:1
   x(k) = ( b(k) -U(k,(k+1):n)*x((k+1):n))/U(k,k);
end

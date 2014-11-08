function   x = tris(U,b)
% Paola Mateos 118798
% Gabriel S�nchez 120141
% Alejandro �lvarez 113638
% Se resuelve el sistema triangular inferior Ux = b usando sustituci�n
% hacia atr�s
% In
% U: Matriz triangular superior.
% b: vector de tama�o n, correspondiente al lado derecho de la ecuaci�n del
% sistema.
%
% Out
% x: vector soluci�n de tama�o n.

n = length(b);
x = zeros(n,1);

x(n) = b(n)/U(n,n);
for k = n-1:-1:1
   x(k) = ( b(k) -U(k,(k+1):n)*x((k+1):n))/U(k,k);
end

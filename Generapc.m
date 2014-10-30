function [G,A,c,b] = Generapc(n,m,tau)
% Genera problemas cuadr?ticos 'aleatorios' de la forma
%   Min  (1/2)*x'*G*x + c'*x
%   s. a.   A*x = b;

% In
% n.- n?mero natural.
% m.- n?mero natural tal que m < n
% tau.- n?mero real positivo menor a 1.
%
% Out
% G.- matriz sim?trica 'aleatoria' de orden n.
%     Se tiene que | G(k,j) | > tau.
% A.- matriz 'aleatoria' de mxn.
%     Se tiene que | A(k,j) | > tau.
% b.-  vector columna de longitud m con puros unos.
% c.- vector columna 'aleatorio' de  longitud.

% generando matrices y vectores
G = triu(rand(n),m); 
A = rand(m,n);
c = rand(n,1);
b = ones(m,1);

% matriz G
for k =1:n
    for j =k:n
        if(abs(G(k,j))<= tau)
            G(k,j) = 0.0;
        end
    end
end
G = (G'*G) + eye(n);

% matriz A
for k =1:m
    for j =1:n
        if(abs(A(k,j))<= tau)
            A(k,j) = 0.0;
        end
    end
end

            
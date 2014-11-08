
function [x] = factible(A,F,b,d)
% esta funcion encuentra un punto factible del conjunto convexo
%        A*x = b
%        F*x >= d,
% por medio de programacion lineal
%In
% A.- matriz de oden mxn
% F.- matriz de orden rxr
% b.- vector columna de orden m.
% d .- vector columna de orden r.
%
%Out
% x.- vector columna de orden, punto factible 
% del conjunto.

% Se plantea el problema de PL
%
% Min       zeros(n,1)'* x + ones(m,1)'* z
% s. a.    A*x = b
%          F*x + z >= d
%             z >= 0
%
% Se resuelve con la funcion en MATLAB:  linprog.m

[m,n] = size(A);
[r,n] = size(F);

 x = zeros(n,1);
 
 %--------------------------------------
 % problema de programacion lineal
 
 
 % vector de la funcion objetivo
 fc = [zeros(n,1); ones(r,1)];
 
% matriz de restricciones de igualdad
 Aeq = [ A  zeros(m,r)];
 
 % matriz de restriccion de desigualdad
 Aiq = [ -F   -eye(r)];
 
 % cotas inferiores y superiores de la variable
 for k = 1:n
    LB(k) = -inf;
    UB(k) = inf;
 end
 
 LB =[LB' ;  zeros(r,1)]';
 
 for k = (n+1):(n+r)
    UB(k) = inf;
 end
 
% Construye un punto factible
R =chol(A*A');
u = R'\b;
v = R\u;
x = A'*v;

z =zeros(r,1);
for j =1:r
    aux = d(j) -F(j,:)*x;
    if (aux > 0)
        z(j) = aux;
    end
    aux = 0;
end

w0 =[x;z];

 % solucion en Matlab
 [w,fw] = linprog(fc,Aiq,-d,Aeq,b,LB,UB,w0);
  
 disp('Valor de la funci?n lineal')
 fw
 

  if (abs(fw) <= 1e-06)
    x = w(1:n);
 else
    disp('El conjunto es numericamente vacio');
    x =[];
 end
 
 
 
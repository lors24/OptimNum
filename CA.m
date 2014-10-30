% Introducir los valores iniciales
%
% Queremos minimizar el problema
% 1/2 * x'*Q*x + c'*x
% s.a. A*x = b
%      F*x >= d
% In:
% Q.- matriz de nxn (sdp)
% c.- vector de nx1
% A.- matriz de mxn (rango(A) = m, m<=n)
% b.- vector de mx1
% F.- matriz de rxn 
% d.- vector de rx1



% Encontrar el conjunto de restricciones activas en xk n 

CA = find(F*xk == d); 
I = [];

m  = length(A); % restricciones de igualdad
R = length(F);  % restricciones de desigualdad
r  = length(I); % restricciones activas de desigualdad

% Resolver el subproblema cuadr?tico

% Si usamos pcnulo encontramos el valor de p, pero faltar?a para resolver
% para lambda
% numero de restricciones en el working set
 
Ak = [A; F(I,:)]; % matriz de restricciones activas
g  = Q*xk + c; % evaluamos el gradiente en xk
b1 = zeros(m+r);

p  = pcnulo( Q, Ak, g, b1);

if (p == 0)
    % calculamos lambda
    
    % encontrar A'
        
    if (lambda <= 0)
        xmin = xk;
    else
        j = find(lambda > 0);
        I = I(I~=j); % quitamos la restriccion j
    end
    
else
    ja = CA(~ismember(CA,I)); %restricciones activas que no estan en el working set
    j = find(F(ja,:)*p <= 0); 
    
    if (isempty(j))  
        jna = find(~ismember(1:R,CA)); % restricciones no activas
        jna1 = find(F(jna,:)*p <= 0);
        [tj,j] = min((d(jna1)-F(jna1,:)*xk) ./ (F(jna1,:)*p));
    
        if (tj < 1)
            xk = xk + tj * p;
            I = [I j]; % agregamos la restriccion que se volvio activa
        else
            xk = xk + p;
        end 
    else
        I = [I f]; %agregamos la restriccion
    end 
end
    

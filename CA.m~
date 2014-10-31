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


%PREPROCESAMIENTO: LA MATRIZ ES REALMENTE LINEALMENTE IDEPENDIENTE

% Encontrar el conjunto de restricciones activas en xk n 
CA = find(F*xk == d); 

% Definir el working set  I c CA
I = [];

m  = length(A); % restricciones de igualdad
R = length(F);  % restricciones de desigualdad
r  = length(I); % restricciones activas de desigualdad EN EL WORKING SET

%RESOLVER SUBPROBLEMA CUADRATICO
%min (1/2)p'Qp + g'p sa Ak*p = 0
 
Ak = [A; F(I,:)]; % matriz de restricciones activas
g  = Q*xk + c; % evaluamos el gradiente en xk
b1 = zeros(m+r);

if (isempty(Ak))
    % Problema sin restricciones
    p = -Q\g
else 
    p  = pcnulo( Q, Ak, g, b1);
end


if (p == 0)
    % CALCULAMOS LAMBDA, UNA VEZ QUE SE TIENE EL xk PODEMOS DESPEJARLA DE
    % LAS CONDICIONES DE PRIMER ORDEN
    lambda = -(Ak*Ak')\(Q*xk + c);
    % encontrar A'
        
    if (lambda <= 0)
        xmin = xk;
    else
        %j = find(lambda > 0); NOTA: ESTO SEG?N YO GENERA UN VECTOR CON
        %TODAS LAS LAMBDAS MAYORES A CERO Y QUEREMOS QUITAR SOLO UNA
        j = find(lambda == max(lambda)); %ENCUENTRA EL M?S NEGATIVO SE PODR?A TENER OTRAS CONDICIONES
        I = I(I~=j); % quitamos la restriccion j
    end
    
else  % p ~= 0 
    ja = find(~ismember(1:R,I)); % RESTRICCIONES FUERA DEL WORKING SET
    j = find(F(ja,:)*p < 0); % RESTRICCIONES DE BLOQUEO
    
    if (~isempty(j))  % HAY RESTRICCIONES DE BLOQUEO
        %jna = find(~ismember(1:R,I)); % RESTRICCIONES FUERA DEL WORKING SET
        %jna1 = find(F(jna,:)*p < 0);
        [tj,ind] = min((d(j)-F(j,:)*xk) ./ (F(j,:)*p));
        ind = ind(1); % PUEDE QUE HAYA DOS ALPHAS IGUALES, TOMAMOS EL PRIMERO
        if (tj < 1)
            xk = xk + tj * p; %damos el paso que nos lleva a la restricci?n de bloqueo
            I = [I j]; % agregamos la restriccion de bloqueo
        else %tj es mayor a 1, das el paso completo, no agregas ninguna restricci?n
             xk = xk + p;
        end 
    else
        xk = xk+p;
    end
   
end
    

function [xmin, iter] = ConjA(Q, c, A, b, F, d)

% realiza una iteracion de ConjA
% Definir el working set  I c CA
I = [];
flag = 0;
xmin = 0;
iter = 0;
[m, ~]  = size(A); % restricciones de igualdad
R = length(d); % restricciones de desigualdad
% restricciones activas de desigualdad EN EL WORKING SET

xk  = factible(A,F,b,d);
%RESOLVER SUBPROBLEMA CUADRATICO
%min (1/2)p'Qp + g'p sa Ak*p = 0

while (flag == 0)
    I = sort(I);
    CA = find(F*xk==d);
    
    Ak = [A; F(I,:)]; % matriz de restricciones activas
    rank(Ak);
    size(Ak);
   %pause(3)
    g  = Q*xk + c; % evaluamos el gradiente en xk
    r  = length(I);
    b1 = zeros(m+r,1);
    
    if (isempty(Ak))
        % Problema sin restricciones
        p = -Q\g;
    else
        p  = pcnulo( Q, Ak, g, b1);
    end
    
    if (norm(p)<1.e-12) %le agregamos una tolerancia
        % CALCULAMOS LAMBDA, UNA VEZ QUE SE TIENE EL xk PODEMOS DESPEJARLA DE
        % LAS CONDICIONES DE PRIMER ORDEN
        lambda = -(Ak*Ak')\(Ak*g);
        % encontrar A'
        
        if (lambda(m+1:end) <= 0)
            xmin = xk;
            flag = 1;
        else
            %TODAS LAS LAMBDAS MAYORES A CERO Y QUEREMOS QUITAR SOLO UNA,
            %las mas positiva
            [~, j] = max(lambda(m+1:end));
            length(lambda)
            j
            size(I)
            I(j) = []; % quitamos la restriccion j
        end
        
    else  % p ~= 0
        
        jna = find(~ismember(1:R,I)); % restricciones no activas
        jna1 = jna(F(jna,:)*p < 0);
        [tj,ind] = min((d(jna1)-F(jna1,:)*xk) ./ (F(jna1,:)*p));
        
        if (tj < 1)
            xk = xk + tj * p;
            I = [I jna1(ind)]; % agregamos la restriccion que se volvio activa
        else
            xk = xk + p;
        end              
        
    end
iter = iter +1;    
end

end










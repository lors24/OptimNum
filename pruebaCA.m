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

% Ingresar datos

Q = 2 * eye(3);
F = [eye(3); - eye(3)];
A = [];
b = [];
c = [-4;-5;-3];
d = [0 0 0 -1 -1 -1]';
xk = [0 0.5 0.5]';

%PREPROCESAMIENTO: LA MATRIZ ES REALMENTE LINEALMENTE IDEPENDIENTE

% Encontrar el conjunto de restricciones activas en xk n 


% Definir el working set  I c CA
I = [];

m  = length(A); % restricciones de igualdad
R = length(F); % restricciones de desigualdad
r  = length(I); % restricciones activas de desigualdad EN EL WORKING SET
flag = 0;

%RESOLVER SUBPROBLEMA CUADRATICO
%min (1/2)p'Qp + g'p sa Ak*p = 0

while (flag == 0)
    
    CA = find(F*xk == d); 
    
    Ak = [A; F(I,:)]; % matriz de restricciones activas
    g  = Q*xk + c; % evaluamos el gradiente en xk
    b1 = zeros(m+r);

    if (isempty(Ak))
     % Problema sin restricciones
        p = -Q\g;
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
            flag = 1;
        else
            %j = find(lambda > 0); NOTA: ESTO SEG?N YO GENERA UN VECTOR CON
            %TODAS LAS LAMBDAS MAYORES A CERO Y QUEREMOS QUITAR SOLO UNA
            j = find(lambda == max(lambda)); %ENCUENTRA EL M?S NEGATIVO SE PODR?A TENER OTRAS CONDICIONES
            I = I(I~=j); % quitamos la restriccion j
        end
    
    else  % p ~= 0 
        ja = CA(~ismember(CA,I)); %restricciones activas que no estan en el working set
        j = find(F(ja,:)*p <= 0); 

        if (isempty(j))   % actualizas x

            jna = find(~ismember(1:R,CA)); % restricciones no activas
            jna1 = jna(F(jna,:)*p <= 0);  
            [tj,j] = min((d(jna1)-F(jna1,:)*xk) ./ (F(jna1,:)*p));
            ind = ind(1); % PUEDE QUE HAYA DOS ALPHAS IGUALES, TOMAMOS EL PRIMERO

            if (tj < 1)
                xk = xk + tj * p;
                I = [I jna1(ind)]; % agregamos la restriccion que se volvio activa
            else
                xk = xk + p;
            end 
        
        else
            I = [I f]; %agregamos la restriccion
        end 
    end
    
end
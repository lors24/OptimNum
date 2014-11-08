
% Para corres una iteracion


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
F = eye(3);
A = [1 1 1];
b = [];
c = [0 0 0]';
d = zeros(3,1);
xk = [0.5 0.5 0]';

%PREPROCESAMIENTO: LA MATRIZ ES REALMENTE LINEALMENTE IDEPENDIENTE

% Encontrar el conjunto de restricciones activas en xk n 


% Definir el working set  I c CA
I = [];

[m, n] = size(A); % restricciones de igualdad
R = length(F); % restricciones de desigualdad
r  = length(I); % restricciones activas de desigualdad EN EL WORKING SET
flag = 0;

%RESOLVER SUBPROBLEMA CUADRATICO
%min (1/2)p'Qp + g'p sa Ak*p = 0
    
    CA = find(F*xk == d);
    I = sort(I);
    Ak = [A; F(I,:)]; % matriz de restricciones activas
    g  = Q*xk + c; % evaluamos el gradiente en xk
    b1 = zeros(m+r,1);

    if (isempty(Ak))
     % Problema sin restricciones
        p = -Q\g;
    else 
        p  = pcnulo( Q, Ak, g, b1);
    end
    
       
    if (norm(p)<1.e-12)
         % CALCULAMOS LAMBDA, UNA VEZ QUE SE TIENE EL xk PODEMOS DESPEJARLA DE
        % LAS CONDICIONES DE PRIMER ORDEN
        lambda = -(Ak*Ak')\(Ak*g);
        % encontrar A'
        
        if (lambda <= 0)
            xmin = xk;
            flag = 1;
        else
            %TODAS LAS LAMBDAS MAYORES A CERO Y QUEREMOS QUITAR SOLO UNA,
            %las mas positiva
            [~, j] = max(lambda); 
            I(j) = []; % quitamos la restriccion j
        end
    
    else  % p ~= 0 
        ja = CA(~ismember(CA,I)); %restricciones activas que no estan en el working set
        j = ja(find(F(ja,:)*p <= 0)); 

        if (isempty(j))   % actualizas x

            jna = find(~ismember(1:R,CA)); % restricciones no activas
            jna1 = jna(F(jna,:)*p < 0);
            [tj,ind] = min((d(jna1)-F(jna1,:)*xk) ./ (F(jna1,:)*p));

            if (tj < 1)
                xk = xk + tj * p;
                I = [I jna1(ind)]; % agregamos la restriccion que se volvio activa
            else
                xk = xk + p;
            end 
            
        else
            I = [I j]; %agregamos la restriccion
        end 
    end
    


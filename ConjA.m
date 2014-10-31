function [I, xk, p, xmin] = ConjA(Q, c, A, b, F, d, xk)

% realiza una iteracion de ConjA
 % Definir el working set  I c CA
I = [];
flag = 0;
xmin = 0;
m  = length(A); % restricciones de igualdad
R = length(F); % restricciones de desigualdad
r  = length(I); % restricciones activas de desigualdad EN EL WORKING SET

%RESOLVER SUBPROBLEMA CUADRATICO
%min (1/2)p'Qp + g'p sa Ak*p = 0
 
while (flag == 0)
    I = sort(I);
    CA = find(F*xk == d); 
    
    Ak = [A; F(I,:)]; % matriz de restricciones activas
    g  = Q*xk + c; % evaluamos el gradiente en xk
    b1 = zeros(m+r,1);

    if (isempty(Ak))
     % Problema sin restricciones
        p = -Q\g;
    else 
        p  = pcnulo( Q, Ak, g, b1);
    end

    if (p == 0)
         % CALCULAMOS LAMBDA, UNA VEZ QUE SE TIENE EL xk PODEMOS DESPEJARLA DE
        % LAS CONDICIONES DE PRIMER ORDEN
        lambda = -(Ak*Ak')\(Ak*g);
        % encontrar A'
        
        if (lambda <= 0)
            xmin = xk;
            flag = 1;
        else
            %j = find(lambda > 0); NOTA: ESTO SEGUN YO GENERA UN VECTOR CON
            %TODAS LAS LAMBDAS MAYORES A CERO Y QUEREMOS QUITAR SOLO UNA
            j = find(lambda == max(lambda)); 
            I = I(I~=j); % quitamos la restriccion j
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
            I = [I j(1)]; %agregamos la restriccion
        end    
        
    end
end

end



    
    



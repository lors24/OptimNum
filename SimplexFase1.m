function [ vinicial ] = SimplexFase1( c, A, b )

[m,n] = size(A);

%PREPARACI�N DEL SUBPROBLEMA
%MIN e'z sa. Ax + Iz = b
%el punto [0;b] es factible
A = [A eye(m)];
x0 = [zeros(n,1); b];
c = [zeros(n,1);ones(m,1)];

%%%%%%%%%%   M�TODO SIMPLEX    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IDENTIFICACI�N DE COLUMNAS B�SICAS Y NO B�SICAS, ASUMIMOS QUE A ES DE
%RANGO COMPLETO Y QUE X0 ES NO DEGENERADO
indb = find(x0 ~= 0);                                                            
indnb = find(x0 == 0);                                                    
%INICIALIZACI�N
%SOL 0,1 INDICA SI EL SIMPLEX ENCONTRO UNA SOLUCI�N, ITER: NUMERO DE
%ITERACIONES QUE TOMO ENCONTRAR LA SOLUCI�N
sol = 0; iter = 0; 

%PROCESO ITERATIVO, CRITERIOS DE PARO: ENCONTRAMOS UNA SOLUCI�N O EL
%PROBLEMA ES NO ACOTADO.
while sol == 0 
    iter = iter+1;
    B = A(:,indb); N = A(:,indnb); cb = c(indb); cn = c(indnb);
    xb = B\b; xn = 0; lambda = B'\cb; sn = cn - N'*lambda;
   
    if sn >= 0
        %ENCONTRAMOS UNA SOLUCI�N
        sol = 1;
        break
    end    
    j = find( sn == min(sn));
    d = B\A(:,j); 
    if d<= 0
        %NO HAY SOLUCI�N, EL PROBLEMA DE MINIMIZACI�N ES NO ACOTADO INFERIORMENTE
        sol = 2;
        break
    end
end

end


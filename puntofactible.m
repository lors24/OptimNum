function xk = puntofactible(A, F, b, d)
% encontrar un punto factible.
% DADO EL PROBLEMA MIN (1/2)x'Qx + c'x sa Ax = b Fx >= d
% Si llevamos las restricciones a la forma estándar de un problema de
% programación lineal podemos encontrar un punto factible
% Recibe A, F, b, d realcionados con las restricciones

% Llevando el problema a forma estándar

 
    % si hay restricciones de ambos tipos va a entrar
    [m,n] = size(A);
    [p,n] = size(F);
    % agregando variables de holgura y separando a x en x+ x-
    AA = [A -A zeros(m, p) ; -F F eye(p)];
    bb = [b;-d];
    % agregamos una z auxiliar
    AA = [ AA eye(m+p)];
    c = [zeros(2*n + p,1); ones(m+p,1)];
    % Resolvemos el problema min e'z sa AAx + z = bb
    X = linprog(c,'','',AA, bb, zeros(2*n + p + m + p),'');
    xk = X(1:n) - X((n+1):(2*n));

end
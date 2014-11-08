%DATOS

% problema 10 de la tarea
Q = 2 * eye(3);
F = [eye(3); - eye(3)];
A = [];
b = [];
c = [-4;-5;-3];
d = [0 0 0 -1 -1 -1]';
xk = [0 0.5 0.5]';

[xmin1, niter1] = ConjA2(Q, c, A, b, F, d)

% problema 9 de la tarea

Q = 2 * eye(2);
F = [1 -2; -1 -2; -1 2; 1 0; 0 1];
A = [];
b = [];
c = [-2 -5]';
d = [-2 -6 -2 0 0]';
xk = [2 0]';

[xmin2, niter2] = ConjA2(Q, c, A, b, F, d)

% simplejo

Q =  2*eye(3);
F = eye(3);
A = [1 1 1];
b = [1];
c = [0 0 0]';
d = zeros(3,1);
%xk = [0.5 0.5 0]';
[xmin3, niter3] = ConjA2(Q, c, A, b, F, d)

% GeneraPC truena

[G,A,F,c,b,d] = Generapc1(500,200,150, 0.5);
[xmin4, niter4] = ConjA2(G,c,A,b,F,d);

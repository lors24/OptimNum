%DATOS

Q = 2 * eye(3);
F = [eye(3); - eye(3)];
A = [];
b = [];
c = [-4;-5;-3];
d = [0 0 0 -1 -1 -1]';
xk = [0 0.5 0.5]';

[xmin, niter] = ConjA(Q, c, A, b, F, d, xk)


%DATOS

Q = 2 * eye(3);
F = [eye(3); - eye(3)];
A = [];
b = [];
c = [-4;-5;-3];
d = [0 0 0 -1 -1 -1]';
xk = [0 0.5 0.5]';

ConjA(Q, c, A, b, F, d, xk)

%[I, xk2, p, xmin] = ConjA(Q, c, A, b, F, d, xk1, I)

%[I, xk3, p, xmin] = ConjA(Q, c, A, b, F, d, xk2, I)

%[I, xk4, p, xmin] = ConjA(Q, c, A, b, F, d, xk3, I)
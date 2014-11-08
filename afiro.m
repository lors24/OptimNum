%Script file: Afiro_quadprog

load afiro

n = length(c);
Q = eye(n);
F = eye(n);
A = full(A);
d = zeros(n,1); % d = lb;

% Min (0.5)*x'*Q*x + c'*x
% s.a. A*x = b
%      F*x >= d

[xmin, niter] = ConjA(Q, c, A, b, F, d)

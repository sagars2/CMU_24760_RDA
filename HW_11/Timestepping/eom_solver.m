function [q_next,dq_next,lambda] = eom_solver(q,dq,h)

m = 1;
g = 9.81;

% q_next = 0;
% dq_next = 0;
% lambda = 0;

M = [m 0;0 m];
N = [0;m*g];
A = [-1 2;1 2];


syms q_next [2,1] real
syms dq_next [2,1] real
syms lambda [2,1] real

a_j = [2*q_next(2)-q_next(1);2*q_next(2)+q_next(1)];
e1 = M*(dq_next-dq)+h*N+A'*lambda == 0;
e2 = q_next-q == h*dq_next;
e3 = a_j >= 0;
e4 = -lambda >= 0;
e5 = a_j.*-lambda == 0;

eqns = [e1,e2,e3,e4,e5];
vars = [q_next,dq_next,lambda];
S = solve(eqns,vars);
q_next = double([S.q_next1; S.q_next2]);
dq_next = double([S.dq_next1; S.dq_next2]);
lambda = double([S.lambda1; S.lambda2]);
end
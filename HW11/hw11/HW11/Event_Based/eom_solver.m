function [q_next,dq_next,lambda] = eom_solver(q,dq,h)

m = 1;
g = 9.81;

M = [m 0;0 m];
N = [0;-m*g];
A = [-1 2;1 2];
a_j = compute_a(in1);

syms q_next[2,1] real
syms dq_next[2,1] real
syms lambda[2,1] real
e1 = M*(dq_next-dq)+h*N+A'*lambda == 0;
e2 = q_next-q == h*dq_next;
e3 = a_j >= 0;
e4 = -lambda >= 0;
e5 = a_j*-lambda == 0;

eqns = [e1,e2,e3,e4,e5,e6];
vars = [q_next,dq_next,lambda];
S = solve(eqns,vars);
end
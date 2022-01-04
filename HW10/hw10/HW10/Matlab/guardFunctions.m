function [constraintFcns,isterminal,direction] = guardFunctions(t,x,contactMode)
m = 1;
g = 9.8;
M = [m 0;0 m];
N = [0;m*g];

dx = x(3);
dy = x(4);
q_dot = [dx;dy];

[a0, a1, a2, a3] = lita_mat(x);
[A0, A1, A2, A3,A13] = A_mat(x);
[dA0, dA1, dA2, dA3, dA13] = dA_mat(x);

if isequal(contactMode,[])
    constraintFcns = [a1;a2;a3];
    isterminal = [1;1;1];
    direction = [-1;-1;-1];
elseif isequal(contactMode,[1])
    constraintFcns = [a2;a3];
    isterminal = [1;1];
    direction = [-1;-1];
elseif isequal(contactMode,[2])
    constraintFcns = [a1];
    isterminal = [1];
    direction = [-1];

elseif isequal(contactMode,[3])
    %find lambda
    isterminal = [1;1;1];
    direction = [-1;-1;1];
    block = inv([M A3.';A3 zeros(1,1)]);
    qddl = block*[-N;-dA3*q_dot];
    lambda = [qddl(3)];
    constraintFcns = [a1;a2;lambda];
elseif isequal(contactMode,[1,3])||isequal(contactMode,[3,1])
    constraintFcns = [a1;a3];
    isterminal = [1;1];
    direction = [-1;-1];

end

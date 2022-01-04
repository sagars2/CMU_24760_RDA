function [dynamicsreturn] = Dynamics(t,x,contactMode)

m = 1;
g = 9.8;
M = [m 0;0 m];
N = [0;m*g];

dx = x(3);
dy = x(4);
q_dot = [dx;dy];

[dA0,dA1,dA2,dA3,dA13] = dA_mat(x);
[A0,A1,A2,A3,A13] = A_mat(x);
lambda = [];
if isequal(contactMode,[])
    ddq = inv(M)*(-N);
    ddx = ddq(1);
    ddy = ddq(2);
    lambda = 0;

elseif isequal(contactMode,[1])
    block = inv([M A1';A1 zeros(1,1)]);
    qddl = block*[-N;-dA1*q_dot];
    ddx = qddl(1);
    ddy = qddl(2);
    lambda = [qddl(3)];
    ddq = [ddx;ddy];
elseif isequal(contactMode,[2])
    block = inv([M A2.';A2 zeros(1,1)]);
    qddl = block*[-N;-dA2*q_dot];
    ddx = qddl(1);
    ddy = qddl(2);
    ddq = [ddx;ddy];
    lambda = [qddl(3)];
elseif isequal(contactMode,[3])
    block = inv([M A3.';A3 zeros(1,1)]);

    qddl = block*[-N;-dA3*q_dot];
    ddx = qddl(1);
    ddy = qddl(2);
    ddq = [ddx;ddy];
    lambda = [qddl(3)];


elseif isequal(contactMode,[1,3]) || isequal(contactMode,[3,1])
    block = inv([M A13.';A13 zeros(2)]);
    qddl = block*[-N;-dA13*q_dot];
    ddx = qddl(1);
    ddy = qddl(2);
    ddq = [ddx;ddy];
    lambda = [qddl(3)];

end
dynamicsreturn = [dx;dy;ddx;ddy];

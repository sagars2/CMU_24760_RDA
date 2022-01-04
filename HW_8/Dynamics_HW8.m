%% Initialization

clear; close; clc;

% States
syms x y theta real 

% Physics
syms w h m real

% Local coordinates
q = [x; y; theta;];

% Constraints
a = [
    y - h*cos(theta) - w*sin(theta);
    y - h*cos(theta) + w*sin(theta);
];

%% Problem 1.1
% We do not need to consider the N, Upsilon or C since N comes from
% potential energy and gravity is always acting. C also comes from mass
% forces which atr always acting as well. None of these parameters affect
% the impact forces and as a result are not considered in this scenario.

% TODO: Compute A matrix
A = sym(zeros(2, 3));
A = [jacobian(a,q)];
Iz = (1/3)*(m*(w^2+h^2));

% TODO: Compute M matrix
M = sym(zeros(3, 3));
M = [m 0 0;0 m 0;0 0 Iz]
%% Problem 1.2

% TODO: Compute dq_plus for the wide block
variables = [m w h x y theta];
values = [1 2 1 2 1 0];
M = subs(M,variables,values)
A = subs(A,variables,values)

dq_minus_wide = [1;-2;-1];
dq_plus_wide = sym(zeros(3, 1));

Big_M = simplify([M A';A zeros(2,2)]);
Big_G = [M*dq_minus_wide;zeros(2,1)];

ans = simplify(inv(Big_M)*Big_G);
dq_plus_wide = simplify([ans(1);ans(2);ans(3)])

P_hat_wide = sym(zeros(2,1));
P_hat_wide = simplify([ans(4);ans(5)])

% TODO: Compute P_hat for the wide block

% TODO: Verify post-impact constraint velocities (A*dq_plus >= 0) for the
% wide block
A_dq_plus_wide = sym(zeros(2,1));
A_dq_plus_wide = A*dq_plus_wide
% TODO: Verify impulses are valid (P_hat <= 0) for the wide block (No
% calculations need to be done for this part as P_hat_wide has already
% been calculated)

%it can be seen that P_hat_wide is negative and is hence <=0

%% Problem 1.3

% TODO: Compute dq_plus for the narrow block
syms x y theta real 
% Physics
syms w h m real
var = [m w h x y theta];
val = [1 1 2 1 2 0];
q_1 = [x;y;theta];
dq_minus_narrow = [2;-1;-1]
A_1 = [jacobian(a,q_1)];
q_1 = subs(q_1,var,val);
A_1 = subs(A_1,var,val)
Iz_1 = subs(Iz,var,val)
M_1 = [m 0 0;0 m 0;0 0 Iz_1];
M_1 = subs(M_1,var,val)

dq_plus_narrow = sym(zeros(3, 1));
Big_M_1 = simplify([M_1 A_1';A_1 zeros(2,2)]);
Big_G_1 = [M_1*dq_minus_narrow;zeros(2,1)];

ans_1 = inv(Big_M_1)*Big_G_1;
dq_plus_narrow = [ans_1(1);ans_1(2);ans_1(3)]
% TODO: Compute P_hat for the narrow block
P_hat_narrow = sym(zeros(2, 1));
P_hat_narrow = [ans_1(4);ans_1(5)]
% TODO: Verify post-impact constraint velocities (A*dq_plus >= 0) are
% either valid or invalid for the narrow block
A_dq_plus_narrow = sym(zeros(2,1));
A_dq_plus_narrow = A_1*dq_plus_narrow
%post-impact constraint velocities constraint is valid.

% TODO: Verify impulses (P_hat <= 0) are either valid or invalid for the
% narrow block (No calculations need to be done for this part as
% P_hat_narrow has already been calculated)

%P_hat is invalid because the top value is greater than 0 
%% Problem 1.4
%setting C1 = 0 since this constraint is non existent and to maintain
%dimensionality
syms x y theta real 
% Physics
syms w h m real
a_2 = [
    y - h*cos(theta) + w*sin(theta);
];

q_2 = [x;y;theta];
dq_minus_narrow = [2;-1;-1]
A_2 = [jacobian(a_2,q_2)];
q_2 = subs(q_2,var,val);
A_2 = subs(A_2,var,val);
Iz_2 = subs(Iz,var,val);
M_2 = [m 0 0;0 m 0;0 0 Iz_2];
M_2 = subs(M_2,var,val);
% TODO: Compute dq_plus for the narrow block with the correct contact mode
dq_plus_correct = sym(zeros(3, 1));
Big_M_2 = simplify([M_2 A_2';A_2 zeros(1,1)]);
Big_G_2 = simplify([M_2*dq_minus_narrow;zeros(1,1)]);
ans_2 = inv(Big_M_2)*Big_G_2;
dq_plus_correct = [ans_2(1);ans_2(2);ans_2(3)]
% TODO: Compute P_hat for the narrow block with the correct contact mode
P_hat_correct = sym(zeros(1, 1));
P_hat_correct = [ans_2(4)]
% TODO: Verify post-impact constraint velocities (A*dq_plus >= 0) are
% valid for the narrow block with the correct contact mode
A_dq_plus_correct = sym(zeros(2,1));
A_dq_plus_correct = A_2*dq_plus_correct
% TODO: Verify impulses (P_hat <= 0) are valid for the narrow block with
% the correct contact mode (No calculations need to be done for this part
% as P_hat_correct has already been calculated)

%P_hat_correct is less than 0. Therefore it is valid
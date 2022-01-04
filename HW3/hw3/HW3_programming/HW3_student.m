format short
clear all; clc; close all;

%% Problem 2.1
% TODO: Writing out by hand the rigid body transformation from the 
% stationary frame to the tool frame in the initial configuration

%% Problem 2.2

% Define symbolic variables
syms q1 q2 q3 q4 q5 q6 real;
q = [q1 q2 q3 q4 q5 q6]';

% Define known frame offsets
lx2 = 320; lx5 = 887; lx6 = 200; lz2 = 680; lz3 = 975; lz4 =200;

% TODO: Define rigid body transformations between successive links
gs1 = [cos(q1) -sin(q1) 0 0; sin(q1) cos(q1) 0 0; 0 0 1 lz2; 0 0 0 1];
g12 = [cos(q2) 0 sin(q2) lx2; 0 1 0 0; -sin(q2) 0 cos(q2) 0;0 0 0 1];
g23 = [cos(q3) 0 sin(q3) 0; 0 1 0 0; -sin(q3) 0 cos(q3) lz3; 0 0 0 1];
g34 = [1 0 0 0; 0 cos(q4) -sin(q4) 0; 0 sin(q4) cos(q4) lz4; 0 0 0 1];
g45 = [cos(q5) 0 sin(q5) lx5; 0 1 0 0; -sin(q5) 0 cos(q5) 0; 0 0 0 1];
g5t = [1 0 0 lx6; 0 cos(q6) -sin(q6) 0; 0 sin(q6) cos(q6) 0; 0 0 0 1];

% TODO: Compute forward kinematics
gst0 = [1 0 0 1407; 0 1 0 0; 0 0 1 1855; 0 0 0 1]

gst = (gs1*g12*g23*g34*g45*g5t)

% TODO: Compute gst(0) with assigning all symbolic variables as 0
gst0_sym = subs(gst, q, [0 0 0 0 0 0]')
% TODO: Compare gst0 to gst0_sym
disp('Element-wise difference between gst0 and gst0_sym: ')
gst0-gst0_sym

%% Problem 2.3

% TODO: Define joint twists in initial configuration
w1 = [0 0 1]';
w2 = [0 1 0]';
w3 = [0 1 0]';
w4 = [1 0 0]';
w5 = [0 1 0]';
w6 = [1 0 0]';

p1 = [0 0 0]';
p2 = [320 0 680]';
p3 = [320 0 1655]';
p4 = [320 0 1855]';
p5 = [1207 0 1855]';
p6 = [1407 0 1855]';

v1 = -cross(w1,p1);
v2 = -cross(w2,p2);
v3 = -cross(w3,p3);
v4 = -cross(w4,p4);
v5 = -cross(w5,p5);
v6 = -cross(w6,p6);

xi1 = [v1;w1];
xi2 = [v2;w2];
xi3 = [v3;w3];
xi4 = [v4;w4];
xi5 = [v5;w5];
xi6 = [v6;w6];

xi1_hat = twist2rbvel(xi1);
xi2_hat = twist2rbvel(xi2);
xi3_hat = twist2rbvel(xi3);
xi4_hat = twist2rbvel(xi4);
xi5_hat = twist2rbvel(xi5);
xi6_hat = twist2rbvel(xi6);
gst_exp = simplify(expm(xi1_hat*q1)*expm(xi2_hat*q2)*expm(xi3_hat*q3)*expm(xi4_hat*q4)*expm(xi5_hat*q5)*expm(xi6_hat*q6)*gst0)

d = subs(gst_exp, q, [0 0 0 0 0 0]');

% TODO: Compare to 2.2
disp('Element-wise difference between gst and gst_exp:')
disp(simplify(gst-gst_exp))

%% Problem 2.4

% TODO: Define goal position as initial configuration +100mm in +y direction 
gDes = [1 0 0 1407; 0 1 0 100; 0 0 1 1855; 0 0 0 1];

% TODO: Define an optimization cost function
F =sum(sum((gst-gDes).^2));
% sum(sum(abs(gst-gDes), "all"))

% Define an optimization problem
fun = matlabFunction(F, 'var', {q}); 

% Call fminunc to solve Inverse Kinematics
options = optimoptions(@fminunc, 'Display', 'iter');
q_sol = fminunc(fun, zeros(6,1), options)


% TODO: Compute the Forward Kinematics using the IK solution
gAchieved = double(subs(gst,q,q_sol))

% TODO: Compare gAchieved with gDes
disp('Element-wise difference between gst and gst_exp:')
difference = (gAchieved - gDes)
%% Problem 2.5

% TODO: Calculate the inverse of Foward Kinematics transformation
gst_inv = inv(gst);

% TODO: Calculate each portion of the spatial Jacobian and get Js
J1 = sym(rbvel2twist(diff(gst,q1)*gst_inv));
J2 = sym(rbvel2twist(diff(gst,q2)*gst_inv));
J3 = sym(rbvel2twist(diff(gst,q3)*gst_inv));
J4 = sym(rbvel2twist(diff(gst,q4)*gst_inv));
J5 = sym(rbvel2twist(diff(gst,q5)*gst_inv));
J6 = sym(rbvel2twist(diff(gst,q6)*gst_inv));
Js = [J1 J2 J3 J4 J5 J6]'

%% Problem 2.6

% TODO: Calculate adjoint of product of exponential map
e2 = sym(expm(xi1_hat*q1));
e3 = sym(expm(xi2_hat*q2));
e4 = sym(expm(xi3_hat*q3));
e5 = sym(expm(xi4_hat*q4));
e6 = sym(expm(xi5_hat*q5));

xi2_prime = sym(tform2adjoint(e2)*xi2);
xi3_prime = sym(tform2adjoint(e2*e3)*xi3);
xi4_prime = sym(tform2adjoint(e2*e3*e4)*xi4);
xi5_prime = sym(tform2adjoint(e2*e3*e4*e5)*xi5);
xi6_prime = sym(tform2adjoint(e2*e3*e4*e5*e6)*xi6);

% TODO: Calculate and compare spatial Jacobian
Js_exp = sym([xi1 xi2_prime xi3_prime xi4_prime xi5_prime xi6_prime]')

% TODO: Compare to 2.6
disp('Element-wise difference between gst and gst_exp:')
simplify(Js_exp - Js)
%% Problem 2.7

% TODO: Define body twist in initial configuration
Vb = [0 1 0 0 0 0]';

% TODO: Convert to spatial twist through adjoint, at initial configuration
Vs = tform2adjoint(gst0).*Vb

% TODO: Compute rank of spatial jacobian (singular if rank < 6)
Js0 = subs(Js_exp,q, [0 0 0 0 0 0]')
rank_Js0 = rank(Js0);
disp('Rank of Js0 is:')
disp(rank_Js0)

disp('No we cannot move the joints in some velocity q(dot) since it loses a rank and has a singularity and the Matrix is singular') 
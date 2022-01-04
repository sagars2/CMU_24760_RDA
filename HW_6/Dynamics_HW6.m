%% Problem 1 Initialization

clear;
clc;

syms q1 q2 dq1 dq2 ddq1 ddq2 real
syms m g l real
syms F tau real
q_g = [q1;q2];
dq_g = [dq1;dq2];
ddq_g = [ddq1;ddq2];
I_o = 1/12*m*l^2;

%% Problem 1.1

% TODO: Kinematic energy
M = [m 0 0;0 m 0;0 0 I_o];
% For computing the kinetic energy, we just need to do 
Jbsl1 = [0 0;l/2 0;1 0];
Jbsl2 = [0 1;(l/2)+q2 0;1 0];
% T_g = simplify(0.5*(dq_g'*(Jbsl1'*M*Jbsl1+Jbsl2'*M*Jbsl2*dq_g))
M_g = (Jbsl1'*M*Jbsl1)+(Jbsl2'*M*Jbsl2)
T_g = 0.5*dq_g'*M_g*dq_g;

% TODO: Potential energy
V_g = simplify((m*g*(l/2)*sin(q1))+m*g*(q2+(l/2))*sin(q1));

% TODO: Lagrangian in generalized coordinates
L_g = simplify(T_g -V_g)

%% Problem 1.2

% TODO: Applied forces to each generalized coordinate
Y_g = sym(zeros(2, 1));

dl_dq_dot = jacobian(L_g,dq_g)
d_dt_dl_dq_dot = simplify(jacobian(dl_dq_dot',q_g)*dq_g + jacobian(dl_dq_dot',dq_g)*ddq_g)

dl_dq = (jacobian(L_g,q_g))';

Upsilon = simplify(d_dt_dl_dq_dot - dl_dq)
Y_g = [tau, F]'
% TODO: Equations of motion by using the Lagrange equations to differentiate the Lagrangian in generalized coordinates
EOM_g1 = simplify(Upsilon-Y_g)

%% Problem 1.3

% TODO: Inertia tensor
M_g = (Jbsl1'*M*Jbsl1)+(Jbsl2'*M*Jbsl2)

% TODO: Coriolis matrix
C_g  = sym(zeros(2, 2));

for i = 1:2
    for j = 1:2
        for k = 1:2
            C_g(i,j) = C_g(i,j)+ (0.5*(diff(M_g(i,j),q_g(k))+diff(M_g(i,k),q_g(j))-diff(M_g(k,j),q_g(i)))*dq_g(k));
        end
    end
end
C_g
                
% TODO: Nonlinear terms
N_g = simplify(jacobian(V_g,q_g))

% TODO: Equations of motion computed directly
EOM_g2 = simplify((M_g*ddq_g)+(C_g*dq_g)+(N_g)'-Y_g)
% simplify(EOM_g1 - EOM_g2)
%% Problem 2 Initialization

% Declare symbolic variables for maximal coordinates and collect into vectors
syms x1 y1 phi_1 x2 y2 phi_2 dx1 dy1 dphi_1 dx2 dy2 dphi_2 ddx1 ddy1 ddphi_1 ddx2 ddy2 ddphi_2 real
syms lambda_1 lambda_2 lambda_3 lambda_4 real
q_m = [x1; y1; phi_1; x2; y2; phi_2];
dq_m = [dx1; dy1; dphi_1; dx2; dy2; dphi_2];
ddq_m = [ddx1; ddy1; ddphi_1; ddx2; ddy2; ddphi_2];
lambda = [lambda_1; lambda_2; lambda_3; lambda_4];
I_o = 1/12*m*l^2;

%% Problem 2.1

% TODO: Constraint function
a_m = sym(zeros(4, 1));
a_m = [x1-(l/2)*cos(phi_1);y1-(l/2)*sin(phi_1);phi_2-phi_1;x2*sin(phi_1)-y2*cos(phi_1)];
% TODO: A as the differential of a
A_m = sym(zeros(4, 6));
A_m = jacobian(a_m, q_m)
%% Problem 2.2

% TODO: Kinematic energy
T_m = sym(0);
T_m = (0.5*m*(dx1^2 + dy1^2))+(0.5*I_o*(dphi_1)^2)+(0.5*m*(dx2^2 + dy2^2))+(0.5*I_o*(dphi_2)^2) %phi1 = phi2
% TODO: Potential energy
V_m = sym(0);
V_m = simplify((m*g*y1)+(m*g*y2))
% TODO: Lagrangian in maximal coordinates
L_m = sym(0);
L_m = simplify(T_m - V_m)

%% Problem 2.3
% TODO: Applied forces in maximal coordinates
Y_m = sym(zeros(6, 1));
Y_m = [-F*cos(phi_1);-F*sin(phi_1);tau; F*cos(phi_1);F*sin(phi_1);0]
%% Problem 2.4

% TODO: Equations of motion by using the Lagrange equations to differentiate the Lagrangian in maximal coordinates
dell_dq_dot = jacobian(L_m,dq_m);
ddt_dldqdot = simplify(jacobian(dell_dq_dot',q_m)*dq_m + jacobian(dell_dq_dot',dq_m)*ddq_m);
dldq = (jacobian(L_m,q_m))';
Upsilon1 = simplify(ddt_dldqdot - dldq);

EOM_m1 = sym(zeros(6, 1));
EOM_m1 = simplify((Upsilon1-Y_m)+A_m'*lambda)
%% Problem 2.5

% TODO: Inertia tensor
M_m = sym(zeros(6, 6));
M_m = jacobian(jacobian(T_m,dq_m),dq_m)
% TODO: Coriolis matrix
C_m  = sym(zeros(6, 6));

for i = 1:6
    for j = 1:6
        for k = 1:6
            C_m(i,j) = C_m(i,j)+ (0.5*(diff(M_m(i,j),q_m(k))+diff(M_m(i,k),q_m(j))-diff(M_m(k,j),q_m(i)))*dq_m(k));
        end
    end
end
C_m
% TODO: Nonlinear terms
N_m = sym(zeros(6, 1));
N_m = simplify(jacobian(V_m,q_m))'
% TODO: Equations of motion computed directly
EOM_m2 = sym(zeros(6, 1));
EOM_m2 = (M_m*ddq_m)+(C_m*dq_m)+ N_m-Y_m+(transpose(A_m)*lambda);

n = simplify(EOM_m2-EOM_m1)
%% Problem 2.6

% TODO: Constraint forces
lambdaVec = sym(zeros(4, 1));
lambdaVec_a = inv(A_m*inv(M_m)*A_m')*(A_m*inv(M_m)*(Y_m-(C_m*dq_m)-N_m));

g = diff(A_m,phi_1)*dphi_1+(diff(A_m,x2))*dx2+diff(A_m,y2)*dy2+diff(A_m,phi_2)*dphi_2;
lambdaVec = simplify(lambdaVec_a+(g*dq_m))

%% Problem 2.7 (Optional)

% TODO: q_m = h(q_g)
h = sym(zeros(6, 1));
h = [(l/2)*cos(q1);(l/2)*sin(q1);q1;(q2+l/2)*cos(q1);(q2+l/2)*sin(q1);q1];
% TODO: dq_m = H*dq_g
H = sym(zeros(6, 2));
H = jacobian(h,q_g);
H = diff(H,q1)*dq1+diff(H,q2)*dq2;
dq_m = H*dq_g;
% TODO: Reconstruct
Mhat = sym(zeros(2, 2));
Mhat = H'*M_m*H;
Chat = sym(zeros(2, 2));
Chat = (H'*M_m*(diff(H,q1)*dq1+diff(H,q2)*dq2))+(H'*C_m*H);
Nhat = sym(zeros(2, 1));
Nhat = H'*N_m;
Yhat = sym(zeros(2, 1));
Yhat = H'*Upsilon1;

% TODO: Equations of motion
EOM_g3 = sym(zeros(2, 1));
EOM_g3 = simplify((Mhat*ddq_g)+(Chat*dq_g)+Nhat-Yhat)
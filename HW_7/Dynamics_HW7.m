%% Initialization

clear;clc;close all

% States
syms theta1 theta2 o_x o_y o_theta dtheta1 dtheta2 do_x do_y do_theta ddtheta1 ddtheta2 ddo_x ddo_y ddo_theta real

% Physics
syms m_l L I_l w m_o I_o g real

% Inputs
syms tau_1 tau_2 real

% Constraint forces
syms lambda_1 lambda_2
lambda = [lambda_1; lambda_2];

% Finger coordinates
theta = [theta1;theta2];
dtheta = [dtheta1;dtheta2];
ddtheta = [ddtheta1;ddtheta2];

% Body coordinates
x = [o_x; o_y; o_theta];
dx = [do_x; do_y; do_theta];
ddx = [ddo_x; ddo_y; ddo_theta];

% Local coordinates
q = [theta;x];
dq = [dtheta; dx];
ddq = [ddtheta; ddx];

%% Problem 1.1

% TODO: Compute inertia matrix in local coordinates Mbar
Mbar = sym(zeros(5, 5));
M_L = [m_l 0 0;
        0 m_l 0;
        0 0 I_l]
M_O = [m_o 0 0;
        0 m_o 0;
        0 0 I_o];
Jbsl1 = [0 0;
        L/2 0;
        1 0]
Jbsl2 = [L*sin(theta2) 0;
        (L/2)+L*cos(theta2) L/2;
        1 1]
M_f = simplify((transpose(Jbsl1)*M_L*Jbsl1)+(transpose(Jbsl2)*M_L*Jbsl2))
M_o = [m_o 0 0;
        0 m_o 0; 
        0 0 I_o]
Mbar = simplify(blkdiag(M_f,M_o))

% TODO: Compute Coriolis matrix in local coordinates Cbar
Cbar  = sym(zeros(5, 5));
% Cbar = [C_f 0;0 C_o];
for i = 1:5
    for j = 1:5
        for k = 1:5
            Cbar(i,j) = Cbar(i,j)+ (0.5*(diff(Mbar(i,j),q(k))+diff(Mbar(i,k),q(j))-diff(Mbar(k,j),q(i)))*dq(k));
        end
    end
end
Cbar
% TODO: Compute nonlinear terms Nbar
Nbar = sym(zeros(5, 1));
% Nbar = [Nf;No];
V_pot = (m_l*g*sin(theta1)*(L/2))+(m_l*g*((L)*sin(theta1)+(L/2)*sin(theta1+theta2)))+m_o*g*o_y;
Nbar = simplify(transpose(jacobian(V_pot,q)))
% TODO: Compute applied force Y
Y = sym(zeros(5, 1));
Y = [tau_1; tau_2;0; 0; 0];

% TODO: Compute A matrix
A = sym(zeros(2, 5));
Bc1 = [1 0 0;
        0 1 0;
        0 0 1; 
        0 0 0; 
        0 0 0;
        0 0 0]
gs1 = [cos(theta1) -sin(theta1) 0 L/2*cos(theta1); 
        sin(theta1) cos(theta1) 0 L/2*sin(theta1);
        0 0 1 0;
        0 0 0 1];
g12 = [cos(theta2) -sin(theta2) 0 L/2+L/2*cos(theta2);
    sin(theta2) cos(theta2) 0 L/2*sin(theta2);
    0 0 1 0;
    0 0 0 1];
gto = [cos(o_theta) -sin(o_theta) 0 (w/2)*cos(-theta1-theta2);
    sin(o_theta) cos(o_theta) 0 (w/2)*sin(-theta1-theta2);
    0 0 1 0
    0 0 0 1];
g2t = [1 0 0 L/2;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1]
gfc = [cos(-theta1-theta2+o_theta) -sin(-theta1-theta2+o_theta) 0 0;
        sin(-theta1-theta2+o_theta) cos(-theta1-theta2+o_theta) 0 0;
        0 0 1 0;
        0 0 0 1]
gst = simplify(gs1*g12*g2t*gfc)

gsc = [cos(o_theta) -sin(o_theta) 0 gst(1,4)
        sin(o_theta) cos(o_theta) 0 gst(2,4);
        0 0 1 0;
        0 0 0 1]
goc = [1 0 0 0;
        0 1 0 -w/2;
        0 0 1 0;
        0 0 0 1];

Ad_gsc1 = tform2adjoint(gsc)
%from notes(ask Ben)
Jsf1 = [0 L*sin(theta1);
        0 -L*cos(theta1);
        0 0;
        0 0;
        0 0;
        1 1];
J_h = simplify(Bc1'*inv(Ad_gsc1)*Jsf1)
J_h_2D = [J_h(1,1) J_h(1,2);
        J_h(2,1) J_h(2,2)]

Ad_goc = tform2adjoint(goc)
G = transpose(Bc1'*inv(Ad_goc))

G_2D = [G((1:2),(1:2));G(6,(1:2))]
gpo = [cos(o_theta) -sin(o_theta) o_x;
        sin(o_theta) cos(o_theta) o_y;
        0 0 1]

Jb_po = sym(zeros(3,3))
for i = 1:3
    temp1 = (simplify(inv(gpo)*diff(gpo,x(i))));
    temp2 = sym(zeros(3,1));
    temp2(1) = temp1(1,3);
    temp2(2) = temp1(2,3);
    temp2(3) = temp1(2,1);
    Jbpo(:,i) = temp2;

end
Jbpo
A = simplify([-J_h_2D G_2D'*Jbpo])

% TODO: Compute dA matrix
dA = sym(zeros(2, 5));
for i = 1:5
    dA = dA+diff(A,q(i))*dq;
end
dA = simplify(dA)
% TODO: Compute equations of motion
EOM = sym(zeros(5, 1));
EOM = simplify((Mbar*ddq)+(Cbar*dq)+(Nbar)+A'*lambda - Y)
%% Problem 1.2

% TODO: Compute acceleration ddq_massive
%Substituting values
% L = 0.1;%m
% g = -9.81;
% m_l = 1;%kg
% tau_1 = 200;%Nm
% tau_2 = 0;%Nm
% I_l = 8.33e-4;%kgm^2
% w = 0.2;%m
% m_o = 24;%kg
% I_o = 0.16;%kg*m^2
% theta1 = pi/2;%degrees
% theta2 = -pi/2;%degrees
% o_x = 0.1;
% o_y = 0.2;
% o_theta = 0;%degrees
% subs(q)
variable = [L,g,m_l,tau_1,tau_2,I_l,w,m_o,I_o,theta1,theta2,o_x,o_y,o_theta,dtheta1,dtheta2,do_x, do_y, do_theta];
value = [0.1,-9.81,1,200,0,8.33e-4,0.2,24,0.16,pi/2,-pi/2,0.1,0.2,0,0,0,0,0,0];

ddq_massive = sym(zeros(5, 1));
% TODO: Compute numerical acceleration ddq_eval_massive
ddq_eval_massive = zeros(5, 1);
b_m_i = inv(simplify([Mbar transpose(A);A zeros(2)]));
b1 = [(Y-(Cbar*dq)-Nbar);-dA*dq];
ddq_massive = (simplify(b_m_i*b1));
ddq_massive_sub = double(subs(ddq_massive,variable,value));
ddq_eval_massive = [ddq_massive_sub(1);ddq_massive_sub(2);ddq_massive_sub(3);ddq_massive_sub(4);ddq_massive_sub(5)];

%% Problem 1.3
% TODO: Compute massless EOM EOM_massless
EOM_massless = sym(zeros(5, 1));
EOM_massless = simplify(subs(EOM,[m_l, I_l],[0,0]))
%% Problem 1.4

% TODO: Compute numerical acceleration ddq_eval_massless
ddq_eval_massless = zeros(5, 1);
variable = [L,g,m_l,tau_1,tau_2,I_l,w,m_o,I_o,theta1,theta2,o_x,o_y,o_theta,dtheta1,dtheta2,do_x, do_y, do_theta];
value = [0.1,-9.81,0,200,0,0,0.2,24,0.16,pi/2,-pi/2,0.1,0.2,0,0,0,0,0,0];

ddq_eval_massless_sub = double(simplify(subs(ddq_massive,variable,value)))

ddq_eval_massless = [ddq_eval_massless_sub(1);ddq_eval_massless_sub(2);ddq_eval_massless_sub(3);ddq_eval_massless_sub(4);ddq_eval_massless_sub(5)]

% TODO: Compute error percentage error_from_massless
error_from_massless = zeros(5, 1);
error_from_massless = (abs((-ddq_eval_massless + ddq_eval_massive))./ddq_eval_massive)*100
%% Problem 1.5

% TODO: Compute the updated A matrix A_frictionless
A_frictionless = sym(zeros(1, 5));
Bc2 = [0;1;0;0;0;0];
G_frictionless_1 = (Bc2'*inv(Ad_goc))';

G_frictionless = [G_frictionless_1(1:2,1);G_frictionless_1(6,1)];
J_h_frictionless = simplify(Bc2'*inv(Ad_gsc1)*Jsf1);

A_frictionless = [-J_h_frictionless G_frictionless'*Jbpo]
% TODO: Compute the updated dA matrix dA_frictionless
dA_frictionless = sym(zeros(1, 5));
for i = 1:5
    dA_frictionless = dA_frictionless+diff(A_frictionless,q(i))*dq(i);
end
dA_frictionless = simplify(dA_frictionless)
% TODO: Compute the updated EOM EOM_frictionless
EOM_frictionless = sym(zeros(5, 1));
EOM_frictionless = simplify((Mbar*ddq)+(Cbar*dq)+(Nbar)+A_frictionless'*lambda_2 - Y)
%% Problem 1.6

% TODO: Compute massless and frictionless EOM EOM_massless_frictionless
EOM_massless_frictionless = sym(zeros(5, 1));
Mbar_massless_frictionless = subs((blkdiag(M_f,M_o)),[m_l,I_l],[0,0])
Nbar_massless_frictionless = subs(Nbar,[m_l,I_l],[0,0])
% Cbar = [C_f 0;0 C_o];
Cbar_massless_frictionless = zeros(5,5)
for i = 1:5
    for j = 1:5
        for k = 1:5
            Cbar_massless_frictionless(i,j) = Cbar_massless_frictionless(i,j)+ (0.5*(diff(Mbar_massless_frictionless(i,j),q(k))+diff(Mbar_massless_frictionless(i,k),q(j))-diff(Mbar_massless_frictionless(k,j),q(i)))*dq(k));
        end
    end
end

Cbar_massless_frictionless 
EOM_massless_frictionless = simplify((Mbar_massless_frictionless*ddq)+(Cbar_massless_frictionless*dq)+(Nbar_massless_frictionless)+A_frictionless'*lambda_2 - Y)
% TODO: Compute the rank of the block matrix
block_massless_frictionless = [Mbar_massless_frictionless A_frictionless';A_frictionless zeros(1,1)]
rank_block = rank(block_massless_frictionless)
%we cannot find the instantaneous acceleration since it is not full rank
%and is hence not invertible. Only one of the links can be taken as
%massless which would allow us to compute instantaneous acceleration.
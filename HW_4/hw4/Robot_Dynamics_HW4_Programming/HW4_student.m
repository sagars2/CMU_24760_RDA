clear; clc; close all;

syms theta_1 theta_2 l_1 l_2 g m real
%% 1.1
% TODO: compute the forward kinematics gst and g_st
g12 = [cos(theta_1) -sin(theta_1) 0 l_1*cos(theta_1); sin(theta_1) cos(theta_1) 0 l_1*sin(theta_1); 0 0 1 0; 0 0 0 1];
g2t = [cos(theta_2) -sin(theta_2) 0 l_2*cos(theta_2);sin(theta_2) cos(theta_2) 0 l_2*sin(theta_2); 0 0 1 0; 0 0 0 1];
g_st=simplify(g12*g2t)
Adg_transpose = (tform2adjoint(g_st))';
F_s = [0; -m*g;0;0;0;-m*g*(l_2*cos(theta_1 + theta_2) + l_1*cos(theta_1))];
% TODO: compute the body wrench Ft as a function of the configuration
F_t=simplify(Adg_transpose*F_s)

%% 1.2
% TODO: compute the body jacobian Jb
Vb_st1 = [l_1*sin(theta_2); l_2+l_1*cos(theta_2); 0; 0;0;1]
Vb_st2 = [0;l_2;0;0;0;1];
Jb = [Vb_st1 Vb_st2]
% gst_inv = inv(g_st);
% J1 = sym(rbvel2twist(diff(g_st,theta_1)*gst_inv));
% J2 = sym(rbvel2twist(diff(g_st,theta_2)*gst_inv));
% Js = [J1 J2]'

% TODO: compute the joint torque tau_pm that counteracts this body wrench
tau_pm = simplify(Jb'*F_t)

%% 1.3
syms m_1 m_2 real

% TODO: compute the joint torques tau_lm that needs to be applied to counteract just the weight of the links themselves
Fs = [0;-m_1*g-m_2*g;0;0;0;-(m_1*g*(l_1/2)*cos(theta_1))-(m_2*g*((l_2/2)*cos(theta_1+theta_2)+l_1*cos(theta_1)))];
Ft = Adg_transpose*Fs;
tau_lm= simplify(Jb'*Ft)

%% 1.4
% TODO: compute the actuator effort 
energy=sym(zeros(1, 1));

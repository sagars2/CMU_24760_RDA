%% 2.6
clear
clc
syms theta_1 theta_2 l1 l2 phi lw ls theta_3 theta_4 a b c d e f g h theta_dot

Bc1 = [1 0 0 0;0 1 0 0;0 0 0 0; 0 0 0 0;0 0 0 0; 0 0 0 0];
Js_s1c1 = [0 l1*sin(theta_1);0 -l1*cos(theta_1);0 0; 0 0;0 0;1 1];
Js_s2c2 = [0 l1*sin(theta_3); 0 -l1*cos(theta_3); 0 0; 0 0; 0 0; 1 1];

gst1 = [cos(theta_1+theta_2) -sin(theta_1+theta_2) 0 l2*cos(theta_1+theta_2)+l1*cos(theta_1);sin(theta_1+theta_2) cos(theta_1+theta_2) 0 l2*sin(theta_1+theta_2)+l1*sin(theta_1);0 0 1 0;0 0 0 1];
gfc1 = [cos(theta_1+theta_2+pi/2-phi) -sin(theta_1+theta_2+pi/2-phi) 0 0;sin(theta_1+theta_2+pi/2-phi) cos(theta_1+theta_2+pi/2-phi) 0 0;0 0 1 0;0 0 0 1]
gsc1 = simplify(gst1*gfc1);
gsc1_inv = inv(gsc1);
adgsc1_inv = simplify(tform2adjoint(gsc1_inv));
J11 = simplify(Bc1'*(adgsc1_inv*Js_s1c1))

gst2 = [cos(theta_3+theta_4) -sin(theta_3+theta_4) 0 l2*cos(theta_3+theta_4)+l1*cos(theta_3);sin(theta_3+theta_4) cos(theta_3+theta_4) 0 l2*sin(theta_3+theta_4)+l1*sin(theta_3);0 0 1 0;0 0 0 1];
gfc2 = [cos((-pi/2)+(theta_3+theta_4+phi)) -sin((-pi/2)+(theta_3+theta_4+phi)) 0 0;sin((-pi/2)+(theta_3+theta_4+phi)) cos((-pi/2)+(theta_3+theta_4+phi)) 0 0;0 0 1 0; 0 0 0 1]
gsc2 = simplify(gst2*gfc2);
gsc2_inv = inv(gsc2);
adgsc2_inv = simplify(tform2adjoint(gsc2_inv));
J22 = simplify(Bc1'*(adgsc2_inv*Js_s2c2));

Jh =[simplify(Bc1'*(adgsc1_inv*Js_s1c1)) zeros(4,2); zeros(4,2) simplify(Bc1'*(adgsc2_inv*Js_s2c2))]

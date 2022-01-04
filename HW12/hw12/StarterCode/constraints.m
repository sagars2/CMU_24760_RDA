function [c,ceq] = constraints(x)

c = [];
ceq = [];
dt = 0.02;
tspan = 0:dt:1.5;
N = length(tspan);
q = x(:,1:2);
dq = x(:,3:4);
tau = x(:,5:6);

qt0 = [-pi/2,0];
dqt0 = [0,0];
qtf = [pi/2,0];
dqtf = [0,0];

q_0 = q(1,:)-qt0;
dq_0 = dq(1,:)-dqt0;

q_f = q(N,:)-qtf;
dq_f = dq(N,:)-dqtf;

dq_e = [];
q_e = [];
ddq_a = [];
% ddq = [];
for i =1:N
    [M,C,N_1,Y] = computeDynamicMatrices(q(i,:)',dq(i,:)',tau(i,:)');
    ddq = inv(M)*(Y-N_1-C*dq(i,:)');
    ddq_a = [ddq_a;ddq'];
end

for i=1:N-1
    dq_e = dq(i,:) + (dt*(ddq_a(i+1,:)+ddq_a(i,:))/2) - dq(i+1,:);
    q_e = q(i,:) + (dt*(dq(i+1,:)+dq(i,:))/2) - q(i+1,:);
    ceq = [ceq,dq_e,q_e];
end

ceq = [ceq,q_0,dq_0,q_f,dq_f];

%q1.3
for j = 1:N
    c = [c,abs(q(j,:))-(3*pi)/4];
end
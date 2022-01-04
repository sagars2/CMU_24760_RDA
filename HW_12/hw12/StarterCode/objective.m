function s_integ = objective(x)
dt = 0.02;
tspan = 0:dt:1.5;
N = length(tspan);
q = x(:,1:2);
dq = x(:,3:4);
tau = x(:,5:6);

s_integ = 0;

for i = 1:N-1
    s_integ = s_integ+(dt*((tau(i,:)*tau(i,:)') + tau(i+1,:)*tau(i+1,:)')/2);
end
s_integ;
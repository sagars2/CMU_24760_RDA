function[CM] = compIV(x,contactMode)
t1 = 0;
x_2 = x(1);
y_2 = x(2);
dx_2 = x(3);
dy_2 = x(4);
[constraintFcns,isterminal,direction] = guardFunctions(t1,x,contactMode);
dq_1 = [dx_2; dy_2];
[a0,a1,a2,a3] = lita_mat(x);
[A0,A1,A2,A3,A13] = A_mat(x);
[dA0,dA1,dA2,dA3,dA13] = dA_mat(x);

A_in = [];
CM_t = [];
m = 1;
g = 9.8;
M = [m 0;0 m];
N = [0;m*g];

dx = x(3);
dy = x(4);
q_dot = [dx;dy];
tol = 0.01;

% Compute reset map into mode
[dq_plus,p_hat] = resetMap(x,contactMode);
% Generate complementarity conditions
if isequal(contactMode,[])
    
    constraintFcns = [a1;a2;a3];
    if abs(a1) <= tol
        CM_t = 1;
        A_in = A1;
    elseif abs(a2)<= tol
        CM_t = 2;
        A_in = A2;
    elseif abs(a3) <= tol
        CM_t = 3;
        A_in = A3;
    end
elseif isequal(contactMode,[1])
    constraintFcns = [a2;a3];
    if abs(a2) <= tol
        CM_t = 2;
        A_in = A2;
    elseif abs(a3) <= tol
        CM_t = [1,3];
        A_in = A3;
    end

elseif isequal(contactMode,[2])
    constraintFcns = [a1];
    if abs(a1) <= tol
        CM_t = 1;
        A_in = A1;
    end

elseif isequal(contactMode,[3])
    CM_t = [];
    A_in = 0;
    
elseif isequal(contactMode,[1,3]) || isequal(contactMode,[3,1])
    constraintFcns = [a1;a3];
    CM_t = [1,3];
    A_in = A13;

end
[dq_plus,p_hat] = resetMap(x,CM_t);

if ((all(p_hat <= tol)) && all(round(A_in*dq_plus,10) >= 0))
        CM = CM_t;
else
    CM = contactMode;
end
end

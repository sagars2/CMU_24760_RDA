function [dq_plus_return,p_hat] = resetMap(x,contactMode)

x_1 = x(1);
y_1 = x(2);
dx_1 = x(3);
dy_1 = x(4);
dq_1 = [dx_1;dy_1];

m = 1;
g = 9.8;
M = [m 0;0 m];
[A0,A1,A2,A3,A13] = A_mat(x);
block = [];

if isequal(contactMode,[])
    dq_plus = dq_1;
    p_hat = 0;
elseif isequal(contactMode,[1])
    block = inv([M A1';A1 zeros(1,1)]);
    A1_dag_tran = [block(1,3);block(2,3)];
    dq_plus = dq_1-(A1_dag_tran)*A1*dq_1;
    p_hat = transpose(A1_dag_tran)*M*dq_1;
elseif isequal(contactMode,[2])
    block = inv([M A2';A2 zeros(1,1)]);
    A2_dag_tran = [block(1,3);block(2,3)];
    dq_plus = dq_1-(A2_dag_tran)*(A2*dq_1); 
    p_hat =transpose(A2_dag_tran)*M*dq_1;
elseif isequal(contactMode,[3])
    block = inv([M A3';A3 zeros(1,1)]);
    A3_dag_tran = [block(1,3);block(2,3)];
    dq_plus = dq_1-A3_dag_tran*A3*dq_1;
    p_hat = transpose(A3_dag_tran)*M*dq_1;
elseif isequal(contactMode,[1,3]) || isequal(contactMode,[3,1])
    block = inv([M A13';A13 zeros(2)]);
    A13_dag_tran = [block(1,3) block(1,4);block(2,3) block(2,4)];
    dq_plus = dq_1-A13_dag_tran*A13*dq_1;
    p_hat = transpose(A13_dag_tran)*M*dq_1;
    
end
dq_plus_return=dq_plus;
end
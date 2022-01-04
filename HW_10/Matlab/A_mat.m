function[A0,A1,A2,A3,A13] = A_mat(x) 
%velocity constraints
A0 = [];
A1 = [0,1]; 
A2 = [1,1];
A3 = [2*x(1)-4, 2*x(2)-2];
A13 = [A1;A3];
end
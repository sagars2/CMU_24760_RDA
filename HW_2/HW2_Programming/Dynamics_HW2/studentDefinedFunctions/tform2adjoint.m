%% tform2adjoint(g) maps the the rigid body transformation g, 
%% in homogeneous coordinates, to the transformation adjoint matrix, Adg.

function Adg=tform2adjoint(g)
    % TODO: construct the 6x6 transformation adjoint matrix Adg from g
    Adg = zeros(6, 6);
    R_ab= [g(1,1) g(1,2) g(1,3); g(2,1) g(2,2) g(2,3); g(3,1) g(3,2) g(3,3)];
    zero = zeros(3,3);
    p_ab = [g(1,4);g(2,4);g(3,4)];
    p_ab_hat = [0 -p_ab(3) p_ab(2);p_ab(3) 0 -p_ab(1);-p_ab(2) p_ab(1) 0]
    d = p_ab_hat * R_ab
    Adg = [R_ab d; zero R_ab]
end

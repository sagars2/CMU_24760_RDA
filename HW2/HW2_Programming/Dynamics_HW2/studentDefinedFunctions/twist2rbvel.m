%% twist2rbvel(xi) maps the 6-vector twist xi to the 4x4 rigid body velocity matrix in homogeneous coordinates xi_hat

function xi_hat=twist2rbvel(xi)
    % TODO: construct the 4x4 rigid body velocity matrix in homogeneous
    % coordinates xi_hat from xi
    w_hat = [0 -xi(6) xi(5); xi(6) 0 -xi(4); -xi(5) xi(4) 0];
    v_skew = [xi(1); xi(2); xi(3)];
    xi_hat = [w_hat v_skew;0 0 0 0];
end
